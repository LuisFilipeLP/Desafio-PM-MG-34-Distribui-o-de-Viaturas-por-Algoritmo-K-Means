#!/usr/bin/env python3
"""
Sistema de apoio à decisão - Policiamento preventivo (PMMG / 39ª CIA).
Backend: K-Means balanceado (algoritmo_kmeans.py) + regras do Cartão-Programa
(Classes 1/2/3, ponto-base 20/40, alternância de quadrículas) + servidor HTTP.

Executar:  python3 patrulha.py [--k 5] [--porta 8000]   ->  http://localhost:8000
Dependências: numpy, scipy (Leaflet/OSM são carregados no navegador).
"""
import argparse, csv, json, math, os, re, threading
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from urllib.parse import urlparse, parse_qs
import numpy as np
from scipy.optimize import linear_sum_assignment
import algoritmo_kmeans as A            # algoritmo original, reutilizado sem alterações

AQUI = os.path.dirname(os.path.abspath(__file__))
A.N_INIT = 10                           # menos reinícios -> recálculo interativo rápido
# Tipos de viatura (premissa configurável; o cenário real virá da organização)
TIPOS = {"V1": "atendimento", "V2": "atendimento"}   # demais = preventiva
RAIO = {1: 1000, 2: 2000, 3: 3000}      # m (classe 2 sem ponto-base = 2 km; classe 3 = 3 km)


def ler_classes():
    """faixas_horarias_classes.csv -> {hora: (classe, modalidade)} (classe 2: 4 ponto-base + 4 patrulha 2 km)."""
    d = {}
    for r in csv.DictReader(open(os.path.join(AQUI, "faixas_horarias_classes.csv"), encoding="utf-8")):
        h = int(r["faixa"][:2]); pb = r["modalidade"].startswith("ponto-base")
        d[h] = (int(r["classe"]), "base" if pb else ("p2" if r["classe"] == "2" else "p3"))
    return d


class Planejador:
    def __init__(self, k):
        self.K = k
        self.H = A.carregar(os.path.join(AQUI, "hexagonos_200m.geojson"))
        self.X, (self.lon0, self.lat0, self.kx, self.ky) = A.projetar(self.H)
        self.classes = ler_classes()
        self.ids = [h["id"] for h in self.H]
        self.cache, self.lock = {}, threading.Lock()
        self.base = None

    def ll(self, xy):
        return [round(self.lat0 + xy[1] / self.ky, 6), round(self.lon0 + xy[0] / self.kx, 6)]

    # ---------- 1) K-Means por hora (idêntico ao main() do algoritmo original) ----------
    def kmeans_24h(self, ativas):
        Ka, H, X = len(ativas), self.H, self.X
        rng = np.random.default_rng(A.SEED); prev = None; out = []
        for h in range(24):
            peso = np.array([x["peso"][h] for x in H], float)
            nivel = np.array([x["nivel"][h] for x in H])
            prio = np.where(nivel >= A.NIVEL_MIN_PRIORITARIO)[0]
            if len(prio) < Ka: prio = np.where(peso > 0)[0]
            if len(prio) < Ka: prio = np.argsort(-peso)[:Ka]
            lab_p, C = A.kmeans_balanceado(X[prio], peso[prio], Ka, rng)
            if prev is None:
                nome = {int(c): ativas[i] for i, c in enumerate(np.argsort(C[:, 0]))}
            else:
                P = np.array([prev[v] for v in ativas])
                r, c = linear_sum_assignment(((C[:, None, :] - P[None, :, :]) ** 2).sum(2))
                nome = {int(i): ativas[j] for i, j in zip(r, c)}
            lab = ((X[:, None, :] - C[None, :, :]) ** 2).sum(2).argmin(1)   # cobertura total (Voronoi)
            lab[prio] = lab_p
            isp = np.zeros(len(H), bool); isp[prio] = True
            prev = {nome[j]: C[j] for j in range(Ka)}
            out.append(dict(peso=peso, prio=isp, viat=np.array([nome[int(l)] for l in lab]), C=prev))
        return out

    # ---------- 2) Seleção da quadrícula-âncora (regra de alternância) ----------
    def escolher_ancoras(self, km, ativas):
        """Âncora = hexágono de maior peso da área da viatura que NÃO foi selecionado na hora
        anterior (nem na seguinte, circularmente) e ainda não foi escolhido por outra viatura na hora."""
        anc = [None] * 24
        ordem = [0] + list(range(1, 24))
        for h in ordem:
            vet = km[h]; proibidos = set()
            for hv in ((h - 1) % 24, (h + 1) % 24):
                if anc[hv]: proibidos |= set(anc[hv].values())
            anc[h] = {}
            for v in ativas:
                idx = np.where(vet["viat"] == v)[0]
                c = np.array(vet["C"][v])
                dist = ((self.X[idx] - c) ** 2).sum(1)
                # prioridade: prioritários > com ocorrência > demais; depois maior peso; depois mais perto do centroide
                chave = sorted(range(len(idx)), key=lambda i: (-int(vet["prio"][idx[i]]), -int(vet["peso"][idx[i]] > 0),
                                                               -vet["peso"][idx[i]], dist[i]))
                for i in chave:
                    hid = self.ids[idx[i]]
                    if hid not in proibidos and hid not in anc[h].values():
                        anc[h][v] = hid; break
        return anc

    # ---------- 3) Cartão-Programa: atividades minuto a minuto ----------
    def rota(self, centro, cand, X, peso, n=6, sep=300):
        """Waypoints: maiores pesos com afastamento mínimo (cobertura), ordenados por vizinho mais próximo."""
        sel = []
        for i in sorted(cand, key=lambda i: (-peso[i], ((X[i] - centro) ** 2).sum())):
            if all(((X[i] - X[j]) ** 2).sum() > sep ** 2 for j in sel):
                sel.append(i)
            if len(sel) >= n: break
        pts, cur, rest = [], centro, list(sel)
        while rest:
            j = min(rest, key=lambda j: ((X[j] - cur) ** 2).sum()); rest.remove(j); pts.append(j); cur = X[j]
        return pts

    def cartoes(self, km, anc, ativas):
        X, ids = self.X, self.ids
        pos = {i: n for n, i in enumerate(ids)}
        cart = {v: [] for v in ativas}
        for h in range(24):
            classe, modal = self.classes[h]
            vet = km[h]
            for vi, v in enumerate(ativas):
                a = pos[anc[h][v]]; R = RAIO[classe] if modal != "base" else 1000
                if modal == "p3": R = 3000
                if modal == "p2": R = 2000
                proprio = np.where(vet["viat"] == v)[0]
                d = np.sqrt(((X - X[a]) ** 2).sum(1))
                dentro = [i for i in proprio if d[i] <= R] or [i for i in np.where(d <= R)[0]]
                way = self.rota(X[a], [i for i in dentro if i != a], X, vet["peso"], n=6)
                A_ = self.ll(X[a]); vll = [self.ll(X[j]) for j in way]; m2 = len(vll) // 2
                if modal == "base":
                    off = 5 + (vi * 9 + h * 7) % 31        # início do ponto-base livre dentro da faixa (05..35)
                    acts = [("patrulha", 0, off, vll[:m2] + [A_]),       # sai da posição anterior, passa pelos pontos e chega à âncora
                            ("base", off, off + 20, [A_]),               # ponto-base de 20 min na quadrícula em destaque
                            ("patrulha", off + 20, 60, [A_] + vll[m2:])]  # 40 min de patrulhamento no total (antes + depois)
                else:
                    acts = [("patrulha", 0, 60, [A_] + vll)]
                cart[v].append(dict(hora=h, classe=classe, modal=modal, ancora=ids[a], raio=R, area=[ids[i] for i in dentro],
                                    atividades=[dict(tipo=t, ini=i0, fim=i1, hex=ids[a], raio=R, rota=r) for t, i0, i1, r in acts]))
        return cart

    def plano(self, indisp):
        ativas = [f"V{i + 1}" for i in range(self.K) if f"V{i + 1}" not in indisp]
        if not ativas: return dict(erro="Nenhuma viatura disponível")
        chave = tuple(ativas)
        with self.lock:
            if chave in self.cache: return self.cache[chave]
        km = self.kmeans_24h(ativas); anc = self.escolher_ancoras(km, ativas)
        res = dict(ativas=ativas, indisp=sorted(indisp), cartoes=self.cartoes(km, anc, ativas),
                   atrib=[[ativas.index(x) for x in km[h]["viat"]] for h in range(24)],
                   prio=[[int(x) for x in km[h]["prio"]] for h in range(24)],
                   centros=[{v: self.ll(c) for v, c in km[h]["C"].items()} for h in range(24)])
        self.verificar(res, anc)
        with self.lock: self.cache[chave] = res
        return res

    def verificar(self, res, anc):
        """Garantias: nenhuma âncora repetida em horas consecutivas; toda viatura tem área; cobertura 100%."""
        for h in range(24):
            a, b = set(anc[h].values()), set(anc[(h + 1) % 24].values())
            assert not (a & b), f"âncora repetida entre {h}h e {(h + 1) % 24}h"
            assert len(set(res["atrib"][h])) == len(res["ativas"]), "viatura sem área"
        res["verificacao"] = "OK"

    def impacto(self, res):
        """Compara com o cenário base (todas disponíveis): hexágonos reatribuídos por hora e equilíbrio de carga."""
        base = self.plano(set()); ativas = res["ativas"]; out = []
        for h in range(24):
            nb = [base["ativas"][i] for i in base["atrib"][h]]; nn = [ativas[i] for i in res["atrib"][h]]
            mud = [i for i in range(len(nb)) if nb[i] != nn[i]]
            pesos = [sum(self.H[i]["peso"][h] for i, a in enumerate(res["atrib"][h]) if a == j) for j in range(len(ativas))]
            out.append(dict(reatribuidos=mud, pesos=pesos, km2=[round(res["atrib"][h].count(j) * 0.025981, 2) for j in range(len(ativas))]))
        return out

    def dados_base(self):
        if self.base is None:
            cia = json.load(open(os.path.join(AQUI, "39_cia.geojson"), encoding="utf-8"))
            pts = []
            try:
                oc = json.load(open(os.path.join(AQUI, "ocorrencias.txt"), encoding="utf-8"))[0]["dados"]
                for o in oc:
                    m = re.match(r"(\d+)/(\d+)/(\d+) (\d+):(\d+)", o["dataHora"])
                    pts.append([round(o["latLon"]["lat"], 5), round(o["latLon"]["lon"], 5), o["natureza"], int(m.group(4))])
            except FileNotFoundError:
                pass
            g = json.load(open(os.path.join(AQUI, "hexagonos_200m.geojson")))
            self.base = dict(cia=cia, ocorrencias=pts, K=self.K, tipos={f"V{i + 1}": TIPOS.get(f"V{i + 1}", "preventiva") for i in range(self.K)},
                             hex=[dict(id=f["properties"]["hex_id"], c=[[p[1], p[0]] for p in f["geometry"]["coordinates"][0]],
                                       nivel=[f["properties"][f"nivel_h{h:02d}"] for h in range(24)],
                                       peso=[f["properties"][f"peso_h{h:02d}"] for h in range(24)]) for f in g["features"]],
                             classes={h: self.classes[h] for h in range(24)})
        return self.base


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--k", type=int, default=5); ap.add_argument("--porta", type=int, default=8000)
    a = ap.parse_args(); P = Planejador(a.k)
    print("Calculando plano base (todas as viaturas)..."); P.plano(set()); print("Pronto.")

    class Hd(BaseHTTPRequestHandler):
        def _send(self, body, ct="application/json"):
            b = body if isinstance(body, bytes) else json.dumps(body, separators=(",", ":")).encode()
            self.send_response(200); self.send_header("Content-Type", ct + "; charset=utf-8")
            self.send_header("Content-Length", str(len(b))); self.end_headers(); self.wfile.write(b)

        def do_GET(self):
            u = urlparse(self.path)
            if u.path == "/api/base": self._send(P.dados_base())
            elif u.path == "/api/plano":
                ind = {s for s in parse_qs(u.query).get("indisp", [""])[0].split(",") if s}
                r = P.plano(ind)
                if "erro" not in r: r = dict(r, impacto=P.impacto(r))
                self._send(r)
            else: self._send(open(os.path.join(AQUI, "index.html"), "rb").read(), "text/html")

        def log_message(self, *a): pass
    print(f"Abra http://localhost:{a.porta}"); ThreadingHTTPServer(("", a.porta), Hd).serve_forever()


if __name__ == "__main__":
    main()
