#!/usr/bin/env python3
"""
K-Means balanceado por hora para distribuição de viaturas (PMMG / 39ª CIA).

Para cada faixa horária (0..23):
  1. Hexágonos PRIORITÁRIOS = hexágonos de maior soma de pesos na hora
     (nível >= NIVEL_MIN_PRIORITARIO, calculado em hexagonos_200m.geojson).
  2. K-Means ponderado (K = nº de viaturas) sobre os centroides dos hexágonos
     prioritários, com peso = soma de pesos do hexágono na hora.
     - Atribuição com RESTRIÇÃO DE TAMANHO: cada viatura recebe floor/ceil(N/K)
       hexágonos prioritários (resolvido como problema de atribuição, Hungarian).
     - Refinamento por trocas para igualar também a SOMA DE PESOS entre viaturas.
     - Várias inicializações (k-means++ ponderado, semente fixa); vence o menor custo.
  3. COBERTURA TOTAL: todo hexágono do território sem prioridade é atribuído ao
     centroide mais próximo (Voronoi) -> nenhuma área fica sem viatura.
  4. Os rótulos V1..VK são casados com os da hora anterior (menor deslocamento)
     para que cada viatura mude o mínimo possível de uma hora para a outra.

Uso:
  python3 algoritmo_kmeans.py [--k 5] [--hex hexagonos_200m.geojson] [--out .]
         [--indisponiveis V2,V4]   # opcional: remove viaturas (cenário dinâmico)
"""
import argparse, csv, json, math, os
import numpy as np
from scipy.optimize import linear_sum_assignment

NIVEL_MIN_PRIORITARIO = 4   # níveis 4 e 5 = hexágonos de maior soma de pesos
N_INIT = 30
LAMBDA_PESO = 1.0           # importância do equilíbrio de pesos vs. compacidade
SEED = 42


def carregar(path):
    g = json.load(open(path, encoding="utf-8"))
    H = []
    for f in g["features"]:
        p = f["properties"]
        H.append(dict(id=p["hex_id"], lat=p["lat"], lon=p["lon"],
                      coords=f["geometry"]["coordinates"][0],
                      peso=[p[f"peso_h{h:02d}"] for h in range(24)],
                      nivel=[p[f"nivel_h{h:02d}"] for h in range(24)],
                      n=[p[f"n_h{h:02d}"] for h in range(24)]))
    return H


def projetar(H):
    lat0 = np.mean([h["lat"] for h in H]); lon0 = np.mean([h["lon"] for h in H])
    kx = 111320 * math.cos(math.radians(lat0)); ky = 110574
    return np.array([[(h["lon"] - lon0) * kx, (h["lat"] - lat0) * ky] for h in H]), (lon0, lat0, kx, ky)


def kmeanspp(X, w, k, rng):
    c = [X[rng.choice(len(X), p=w / w.sum())]]
    for _ in range(1, k):
        d2 = np.min([((X - ci) ** 2).sum(1) for ci in c], axis=0)
        p = d2 * w
        c.append(X[rng.choice(len(X), p=p / p.sum())] if p.sum() > 0 else X[rng.integers(len(X))])
    return np.array(c)


def atribuicao_balanceada(X, C):
    """Cada centroide recebe floor(N/K) ou ceil(N/K) pontos, minimizando a distância².
    Colunas: q vagas "regulares" por viatura + 1 vaga "extra" por viatura;
    (K - r) linhas fantasmas ocupam vagas extras, sobrando r extras para pontos reais."""
    n, k = len(X), len(C)
    q, r = divmod(n, k)
    slots = [j for j in range(k) for _ in range(q)] + (list(range(k)) if r > 0 else [])   # regulares + extras
    extra0 = k * q
    d2 = ((X[:, None, :] - C[None, :, :]) ** 2).sum(2)
    cost = d2[:, slots]
    fant = k - r if r > 0 else 0
    if fant:
        BIG = 1e12
        dummy = np.full((fant, len(slots)), BIG); dummy[:, extra0:] = 0.0
        cost = np.vstack([cost, dummy])
    rr, cc = linear_sum_assignment(cost)
    lab = np.empty(n, int)
    for i, c in zip(rr, cc):
        if i < n:
            lab[i] = slots[c]
    return lab


def custo(X, w, lab, C, d0):
    k = len(C)
    d = sum(((X[lab == j] - C[j]) ** 2).sum() for j in range(k))
    S = np.array([w[lab == j].sum() for j in range(k)])
    return d / (len(X) * d0) + LAMBDA_PESO * (((S - S.mean()) / S.mean()) ** 2).sum()


def centroides(X, w, lab, k, C_old):
    C = C_old.copy()
    for j in range(k):
        m = lab == j
        if m.any():
            C[j] = (X[m] * w[m, None]).sum(0) / w[m].sum()
    return C


def refinar_trocas(X, w, lab, C, d0):
    """Trocas de pares entre viaturas que reduzem o custo (mantém tamanhos iguais)."""
    melhor = custo(X, w, lab, C, d0)
    melhorou = True
    while melhorou:
        melhorou = False
        for i in range(len(X)):
            for j in range(i + 1, len(X)):
                if lab[i] == lab[j]:
                    continue
                lab[i], lab[j] = lab[j], lab[i]
                c = custo(X, w, lab, C, d0)
                if c < melhor - 1e-12:
                    melhor, melhorou = c, True
                else:
                    lab[i], lab[j] = lab[j], lab[i]
    return lab


def kmeans_balanceado(X, w, k, rng):
    best = None
    for _ in range(N_INIT):
        C = kmeanspp(X, w, k, rng)
        lab = atribuicao_balanceada(X, C)
        d0 = max(1.0, np.mean(((X - C[lab]) ** 2).sum(1)))
        for _ in range(50):
            C = centroides(X, w, lab, k, C)
            lab = atribuicao_balanceada(X, C)
            lab = refinar_trocas(X, w, lab, C, d0)
            C2 = centroides(X, w, lab, k, C)
            if np.allclose(C2, C):
                break
            C = C2
        c = custo(X, w, lab, C, d0)
        if best is None or c < best[0]:
            best = (c, lab.copy(), C.copy())
    return best[1], best[2]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--k", type=int, default=5)
    ap.add_argument("--hex", default="hexagonos_200m.geojson")
    ap.add_argument("--out", default=".")
    ap.add_argument("--indisponiveis", default="", help="ex.: V2,V4")
    a = ap.parse_args()

    H = carregar(a.hex)
    X, (lon0, lat0, kx, ky) = projetar(H)
    K = a.k
    indisp = {s.strip() for s in a.indisponiveis.split(",") if s.strip()}
    ids_viaturas = [f"V{i+1}" for i in range(K)]
    ativas = [v for v in ids_viaturas if v not in indisp]
    Ka = len(ativas)
    if Ka == 0:
        raise SystemExit("Nenhuma viatura disponível.")
    rng = np.random.default_rng(SEED)
    area_hex_km2 = 2.598076 * 100 ** 2 / 1e6

    linhas, resumo, centros_json = [], [], {}
    prev = None  # centroides da hora anterior {viatura: xy}

    for h in range(24):
        peso = np.array([x["peso"][h] for x in H], float)
        nivel = np.array([x["nivel"][h] for x in H])
        prio = np.where(nivel >= NIVEL_MIN_PRIORITARIO)[0]
        if len(prio) < Ka:                       # poucas ocorrências na hora
            prio = np.where(peso > 0)[0]
        if len(prio) < Ka:                       # fallback: top-K por peso
            prio = np.argsort(-peso)[:Ka]

        lab_p, C = kmeans_balanceado(X[prio], peso[prio], Ka, rng)

        # casar clusters com as viaturas da hora anterior (mínimo deslocamento)
        if prev is None:
            ordem = np.argsort(C[:, 0])          # V1..Vk de oeste para leste
            nome = {int(c): ativas[i] for i, c in enumerate(ordem)}
        else:
            P = np.array([prev[v] for v in ativas])
            cost = ((C[:, None, :] - P[None, :, :]) ** 2).sum(2)
            r, c = linear_sum_assignment(cost)
            nome = {int(i): ativas[j] for i, j in zip(r, c)}

        # cobertura total: demais hexágonos -> centroide mais próximo
        d2 = ((X[:, None, :] - C[None, :, :]) ** 2).sum(2)
        lab_all = d2.argmin(1)
        lab_all[prio] = lab_p                     # prioritários mantêm a atribuição balanceada
        is_prio = np.zeros(len(H), bool); is_prio[prio] = True

        # verificações de integridade
        assert len(set(lab_all)) == Ka, "viatura sem área"
        cnt = np.bincount(lab_p, minlength=Ka)
        assert cnt.max() - cnt.min() <= 1, "prioritários desbalanceados"

        prev = {nome[j]: C[j] for j in range(Ka)}
        centros_json[h] = {nome[j]: [round(lat0 + C[j][1] / ky, 6), round(lon0 + C[j][0] / kx, 6)] for j in range(Ka)}

        for i, x in enumerate(H):
            linhas.append([h, x["id"], nome[int(lab_all[i])], int(is_prio[i]),
                           int(peso[i]), int(nivel[i]), x["n"][h]])
        for j in range(Ka):
            m = lab_all == j; mp = m & is_prio
            resumo.append([h, nome[j], int(m.sum()), round(m.sum() * area_hex_km2, 2),
                           int(mp.sum()), int(peso[mp].sum()), int(peso[m].sum()),
                           round(lat0 + C[j][1] / ky, 6), round(lon0 + C[j][0] / kx, 6)])

    os.makedirs(a.out, exist_ok=True)
    with open(os.path.join(a.out, "kmeans_hexagono_viatura_por_hora.csv"), "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["hora", "hex_id", "viatura", "prioritario", "soma_pesos", "nivel", "n_ocorrencias"])
        w.writerows(linhas)
    with open(os.path.join(a.out, "kmeans_resumo_viaturas_por_hora.csv"), "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["hora", "viatura", "hexagonos_area", "area_km2", "hex_prioritarios",
                    "peso_prioritarios", "peso_total_area", "centro_lat", "centro_lon"])
        w.writerows(resumo)
    json.dump(dict(viaturas=ativas, centros=centros_json,
                   atrib=[[l[2] for l in linhas[h * len(H):(h + 1) * len(H)]] for h in range(24)],
                   prio=[[l[3] for l in linhas[h * len(H):(h + 1) * len(H)]] for h in range(24)],
                   hex_ids=[x["id"] for x in H]),
              open(os.path.join(a.out, "kmeans_atribuicao.json"), "w"), separators=(",", ":"))
    print(f"OK: K={Ka} viaturas ({', '.join(ativas)}), {len(H)} hexágonos, 24 horas.")


if __name__ == "__main__":
    main()
