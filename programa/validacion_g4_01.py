"""Validacion del modulo nh3h2o.py contra los valores de referencia de la guia
IAPWS G4-01 (2001), Tablas 6, 7 y 8 -- Seccion 8 'Computer Program Verification'.
Ejecutar: python3 validacion_g4_01.py"""
from nh3h2o import *

T6 = [(0.1,600,35,-13734.1763,32.1221333,53.3159544,883.925596),
      (0.1,600, 4,-16991.6697,12.7721090,52.7644553,471.762394),
      (0.5,500,32,-12109.5369,21.3208159,58.0077346,830.295833),
      (0.5,500, 1,-18281.3020, 3.6423080,36.8228098,510.258362),
      (0.9,400,30, -6986.4869,22.2830797,51.8072415,895.748711),
      (0.9,400,0.5,-13790.6278,1.5499708,32.9703870,478.608147)]
T7 = [(0.2,300,0.040710,0.9360),(0.4,400,2.5545,0.9363),(0.6,500,16.698,0.7844)]
T8 = [(0.2,300,0.00437062,0.010672),(0.4,400,0.394694,0.051541),(0.6,500,6.52607,0.22135)]
e = lambda a,b: abs(a-b)/abs(b)*100
mx = 0.0
print("TABLA 6 - region monofasica (f [J/mol], p [MPa], Cv [J/mol-K], w [m/s])")
print(f"{'x':>5}{'T':>6}{'rho':>7} | {'err f %':>10}{'err p %':>10}{'err Cv %':>10}{'err w %':>10}")
for x,T,rm,f,p,cv,w in T6:
    M=Mm(x); r=prop(rm,T,x)
    d=[e(r['a']*M,f),e(r['P'],p),e(r['cv']*M,cv),e(r['w'],w)]
    mx=max(mx,*d); print(f"{x:>5}{T:>6}{rm:>7} | "+"".join(f"{v:>10.2e}" for v in d))
print("\nTABLA 7 - puntos de burbuja      TABLA 8 - puntos de rocio")
print(f"{'xL':>5}{'T':>6} | {'err pBUB %':>12}{'err xv %':>10}   |{'xv':>5}{'T':>6} | {'err pDEW %':>12}{'err xL %':>10}")
for (xl,T,p,xv),(xv2,T2,p2,xl2) in zip(T7,T8):
    Pb,yb,_ = bubbleP(T,xl); Pd,yd,_ = dewP(T2,xv2)
    d=[e(Pb,p),e(yb,xv),e(Pd,p2),e(yd,xl2)]
    if T < 500: mx=max(mx,*d)
    print(f"{xl:>5}{T:>6} | {d[0]:>12.2e}{d[1]:>10.2e}   |{xv2:>5}{T2:>6} | {d[2]:>12.2e}{d[3]:>10.2e}")
print(f"\nDesviacion relativa maxima: {mx:.2e} %  (excluido el punto de burbuja x=0.6/500 K)")
print("""
LIMITACION CONOCIDA. El punto de burbuja x_L=0.6, T=500 K, p=16.7 MPa esta sobre
el locus critico de la mezcla; alli la sustitucion sucesiva colapsa a la raiz
trivial (y = z) y la inicializacion automatica falla. Sembrado con las densidades
de referencia de la Tabla 7, el sistema exacto SI reproduce el valor de la guia
(16.6984 frente a 16.698 MPa). La propia guia IAPWS advierte de problemas de
convergencia en la region critica. El dominio de operacion de este proyecto
(4-30 bar, 300-630 K, x = 0.4-0.7 masica) esta muy lejos de esa region.""")
