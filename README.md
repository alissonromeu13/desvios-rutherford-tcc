# Scripts das figuras do TCC

Alisson Romeu de Sousa, *Desvios da lei de espalhamento alfa de Rutherford* (IFCE Crateús, 2026).

Cada figura do TCC tem um script, numerado como no documento. Rodar todos reproduz, em `figuras/`, os 18 PDFs com os mesmos nomes e pastas usados no Overleaf. Os scripts de `numeros_do_texto/` imprimem os números citados no texto que não aparecem nas próprias figuras.

## Como rodar

1. Instale Python 3.13 e os pacotes de `requirements.txt` (`pip install -r requirements.txt`).
2. Instale uma distribuição LaTeX com XeLaTeX e as fontes Latin Modern (MiKTeX ou TeX Live). O estilo usa o backend `pgf`, então o texto das figuras é composto pelo XeLaTeX com as mesmas fontes do documento.
3. Rode `python rodar_todas.py` para gerar tudo (cerca de 1 min), ou `python rodar_todas.py fig17` para gerar uma só figura.

Os PDFs saem em `figuras/<pasta>/<nome>.pdf`. Se o `pdftoppm` (Poppler) estiver instalado, sai também um PNG de pré-visualização.

## Módulos comuns

| Arquivo | Conteúdo |
|---|---|
| `estilo.py` | Fontes, cores, vírgula decimal (ABNT) e a função `salvar()` |
| `scatter_core.py` | Constantes (CODATA 2018), funções de Coulomb $F_l$ e $G_l$ (algoritmo de Steed; Barnett et al. 1974), integração de Numerov, matriz $S$, amplitudes, Woods-Saxon, Coulomb de esfera uniforme e conversão laboratório → centro de massa |

## Figuras

| Fig. | Seção | Script | Arquivo no Overleaf (`figuras/…`) | Origem dos parâmetros e dados |
|---|---|---|---|---|
| 1 | 2.1.1 | `fig01_rutherford_trajetoria.py` | `2/fig_rutherford_trajetoria.pdf` | Esquema analítico (órbita de Rutherford) |
| 2 | 2.1.2 | `fig02_rutherford_feixe.py` | `2/fig_rutherford_feixe.pdf` | Esquema analítico |
| 3 | 2.3.1 | `fig03_regimes_sommerfeld.py` | `4.1/fig01_regimes_sommerfeld_v2.pdf` | Constantes CODATA; $R_{\text{toque}} = 1{,}2\,(A_\alpha^{1/3} + A^{1/3})$ fm (Krane, p. 48 e 480); fronteira $\eta = 4$ por convenção (Bohr 1948; Mott e Massey, p. 126) |
| 4 | 2.3.2 | `fig04_serie_born.py` | `2/fig_serie_born.pdf` | Esquema |
| 5 | 2.4.1 | `fig05_esfera_homogenea.py` | `5.1/fig_5_1_esfera_homogenea.pdf` | Analítico, $F = 3j_1(x)/x$ |
| 6 | 2.4.2 | `fig06_exponencial_dipolo.py` | `5.2/fig_5_2_exponencial_dipolo.pdf` | Analítico, $F = (1+x^2)^{-2}$ |
| 7 | 2.4.3 | `fig07_gaussiana.py` | `5.3/fig_5_3_gaussiana.pdf` | Analítico, $F = e^{-x^2/4}$ |
| 8 | 2.4.4 | `fig08_validade_born_gaussiano.py` | `5.3/fig02_validade_born_gaussiano_v2.pdf` | Poço gaussiano ilustrativo ($a = 2$ fm; $V_0 = -2$ e $-20$ MeV; massa do nêutron), escolhido para $\xi \ll 1$ e $\xi \sim 1$ |
| 9 | 2.4.5 | `fig09_fatores_forma_realistas.py` | `5.4/fig10_fatores_forma_realistas_v4.pdf` | 2pF do $^{207}$Pb, $c = 6{,}62$ fm e $d = 0{,}546$ fm (difusividade) (De Vries et al. 1987, p. 510, Tab. I) |
| 10 | 4.1.2 | `fig10_born_bieler.py` | `Potencial_combinado_e_superposicao_de_amplitudes/fig10_born_bieler.pdf` | Bieler (1924), p. 445, Tabs. I–III; peso da Eq. 1, p. 438; energia pela âncora do Apêndice F |
| 11 | 4.2.2 | `fig11_onda_s.py` | `6.3/fig03_onda_s_comprimento_espalhamento_v2.pdf` | Analítico; pontos n–p com $a$ e $r_e$ de Babenko e Petrov (2007), p. 672, Eq. 17 |
| 12 | 4.2.3 | `fig12_ondas_parciais_ressonancia.py` | `6.3/fig04_ondas_parciais_ressonancia_v3.pdf` | Analítico; $K_0R = 5{,}30$ (escolha justificada em `numeros_do_texto/n12_escolha_K0R.py`) |
| 13 | 4.2.4.2 | `fig13_disco_negro.py` | `6.3/fig05_disco_negro_difracao_v3.pdf` | Analítico ($kR = 10$ e 30) |
| 14 | 4.2.4.2 | `fig14_sombra_disco_negro.py` | `6.3/fig12_sombra_disco_negro.pdf` | Analítico (Fresnel e Fraunhofer) |
| 15 | 4.4.1 | `fig15_poco_complexo.py` | `6.4/fig06_poco_complexo_absorcao_v4.pdf` | Parâmetros ilustrativos do tipo nêutron + $^{27}$Al ($\mu = 906$ MeV, $R = 3{,}6$ fm, $V_0 = 40$ MeV), declarados no texto |
| 16 | 4.4.2 | `fig16_woods_saxon_neutro.py` | `6.4/fig07_woods_saxon_neutro_v4.pdf` | Koning e Delaroche (2003): Eqs. 1–3 e 7, Tabs. 5 e 6 (sem spin-órbita) |
| 17 | 4.4.3 | `fig17_alfa_aluminio_modelo_optico.py` | `6.4/fig08_alfa_al_modelo_optico_bieler_v3.pdf` | McFadden e Satchler (1966), p. 198 ($V = 185$, $W = 25$ MeV, $r_0 = 1{,}40$, $a = 0{,}52$, $r_C = 1{,}3$ fm); Bieler (1924), Tabs. I e II; Apêndice F |
| 18 | 4.4.3 | `fig18_alfa_ouro_energia.py` | `6.4/fig09_alfa_au_dependencia_energia.pdf` | Mesmo potencial; Farwell e Wegner (1954), p. 1216, Tab. I; fatores cinemáticos do Apêndice E |

## Números citados no texto

| Script | O que imprime |
|---|---|
| `numeros_do_texto/n03_raio_de_corte.py` | Raio de corte da condição de fase, $2\eta K_0(r_{\text{lim}}/R) = 1$, para $\alpha$ + H, ouro e elétrons de 10 keV no hidrogênio, comparado às formas assintóticas; critério de altas energias $2\eta\ln(kr_B)$ e fase na escala de $\bar\lambda$ (Seção 2.3.3) |
| `numeros_do_texto/n10_fator_forma_carga_bieler.py` | $\lvert F(q)\rvert^2$ da esfera homogênea nos ângulos efetivos dos diafragmas de Bieler, comparado às razões medidas; limite $g < 0{,}94$ para o mínimo do modelo de poço em $qR \approx 2{,}74$; comparação com as distribuições exponencial e gaussiana de mesmo raio quadrático médio (Seção 4.1.2) |
| `numeros_do_texto/n12_escolha_K0R.py` | Varredura de $K_0R$ que levou à ressonância $l = 3$ estreita da Figura 12 |
| `numeros_do_texto/n13_escolha_Delta.py` | Oscilação de σ_reação com o corte abrupto em *l* e com a borda suave Δ = 0,5 da Figura 13(b) |
| `numeros_do_texto/n17_alfa_aluminio.py` | $\sigma_{\text{reação}}$ em 24,7 MeV para os três potenciais do Al (McFadden e Satchler, Tab. 3, p. 187); $\chi^2$ dos seis pontos de Bieler; funções de excitação em 90° e 135° comparadas a Rutherford e Chadwick (1925) e a Farwell e Wegner (1954, p. 1212); teste com absorção puramente de superfície |
| `numeros_do_texto/n18_alfa_ouro.py` | Ausência de desvio no ouro; energias $E_{1/4}$ do modelo e de Farwell e Wegner; inclinação $K$; fatores cinemáticos |

Os scripts `fig10`, `fig17` e `fig18` também imprimem, ao rodar, os demais números das Seções 4.1.2 e 4.4.3.

## Histórico

- **1.1 (22 nov. 2026):** Figura 1 com vetores em negrito, $\varphi'$, $\beta$, $r_{\text{min}}$ e o rótulo "partícula alfa"; Figura 2 com o arco $r\,|\mathrm{d}\theta|$, o segmento $r\,\mathrm{sen}\,\theta$, a faixa $\mathrm{d}A$ sombreada e setas nas trajetórias incidentes; Figura 3 com as partículas alfa de Bieler (7,35 MeV) no ouro e no alumínio, rótulos com fundo branco e o regime "perturbativo coulombiano"; eixos $\mathrm{d}\sigma/\mathrm{d}\Omega$ com d reto (Figuras 8, 13 e 16); legendas das Figuras 8, 9 e 12; novo `n03_raio_de_corte.py`, com o ouro em 7,35 MeV ($\eta \approx 18{,}4$); `n10` compara a esfera com as distribuições exponencial e gaussiana de mesmo raio quadrático médio; números novos impressos por `fig12` e `fig17`.
- **1.0 (16 nov. 2026):** versão citada na primeira entrega.

## Fora do TCC

`extras/argand_ressonancia_nao_usada.py` gera o diagrama de Argand da ressonância $l = 3$. A figura não entrou no texto e foi guardada como possível apêndice.

## Convenções

- Energias em MeV, comprimentos em fm, ângulos de laboratório quando o índice "lab" aparece.
- A energia das partículas alfa de Bieler vem da lei de Geiger, $E \propto R^{2/3}$, com $E(7{,}0\ \text{cm}) = 7{,}644$ MeV (Apêndice F).
- As médias sobre os diafragmas de Bieler usam o peso $\operatorname{cossec}^3(\theta/2)$ da geometria do experimento (Bieler 1924, p. 438, Eq. 1), nas Figuras 10 e 17.
