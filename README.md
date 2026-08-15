# Producción de ¹³¹I por activación neutrónica de TeO₂ — simulaciones

**Autor:** Adrian Navarro Minaya · [adriannavarrominaya@gmail.com](mailto:adriannavarrominaya@gmail.com)

**Tutor:** Oscar Luis Cabellos de Francisco · [oscar.cabellos@upm.es](mailto:oscar.cabellos@upm.es)

**Departamento:** Departamento de Ingeniería Energética (Área Nuclear)

**Centro:** Universidad Politécnica de Madrid

**Fecha:** Agosto 2026

**Códigos de simulación:** ACAB 2008 (UPM — Activation code) · COLLAPS (preprocesador de secciones eficaces)

---

Material de cálculo del Trabajo de Fin de Grado *«Optimización y caracterización
física de la producción de ¹³¹I a partir de la activación con neutrones de óxido de
telurio»*.

Este repositorio contiene **todas las simulaciones ejecutadas**, con sus ficheros de
entrada y de salida sin editar, y los scripts de verificación independiente que se
usaron para contrastarlas. No contiene la suite de aplicaciones, que vive en su
propio repositorio.

**El motor de cálculo es ACAB** (*Activation Abacus Code*, Sanz, Cabellos y
García-Herranz, 2008) y **no se ha modificado**. Lo que aquí se publica son sus
entradas, sus salidas y las comprobaciones hechas por fuera.

---

## Qué se calcula

La captura neutrónica sobre el ¹³⁰Te —el 33,8 % del telurio natural— produce ¹³¹I por
dos caminos:

```
¹³⁰Te ─(n,γ)→ ¹³¹Te  (25 min) ─β⁻→ ¹³¹I  (8,0252 d)
¹³⁰Te ─(n,γ)→ ¹³¹ᵐTe (30 h)   ─β⁻/IT→ …
```

El canal metaestable recibe el 4,2 % de las capturas y es el que gobierna la forma de
la curva: **el máximo de actividad no está en el instante de apagar el reactor**.

---

## Estructura

```
Conjunto Simulaciones 1 - Verificacion/
Conjunto Simulaciones 2 - Simulacion de Referencia/
Conjunto Simulaciones 2 - Barridos/
Conjunto Simulaciones 2 - Análisis de cadenas/
Scripts/
```

### `Conjunto Simulaciones 1 - Verificacion`

Réplica de los cuatro experimentos del trabajo de referencia (Haffner, Miller y
Morris, 2019), cada uno con varias bibliotecas de datos nucleares. Los cuatro
comparten blanco y difieren en el historial temporal y el flujo:

| Experimento | Irradiación | Flujo total [n·cm⁻²·s⁻¹] |
|---|---|---|
| 1 | 10,0008 s | 6,67×10¹³ |
| 2 | 20 min | 6,67×10¹³ |
| 3 | 30 min | 8,51×10¹³ |
| 4 | 150 h | 1,0038×10¹⁴ |

### `Conjunto Simulaciones 2 - Simulacion de Referencia`

La simulación sobre la que se construyen los barridos y el análisis de cadenas.
Punto de trabajo del cuarto experimento, biblioteca de 211 grupos con decaimiento de
2025, y una malla temporal de **79 nodos de irradiación y 120 de enfriamiento**.

Sus 79 nodos de irradiación **son en sí mismos el barrido de historial temporal**:
cada uno da el inventario para su tiempo de irradiación, y desde cualquiera de ellos
se puede propagar analíticamente el enfriamiento.

### `Conjunto Simulaciones 2 - Barridos`

Tres barridos paramétricos, cada uno con su `sweep_manifest.json`:

| Barrido | Puntos | Qué cambia del `inp.5` |
|---|---|---|
| Flujo | 5, de 5×10¹² a 1×10¹⁷ | **`XNORM`, bloque 9** — no el flujo del bloque 3 |
| Masa | 6, de 0,123 a 1000 g | Bloque 5 |
| Espectro | 9 espectros de instalación | **El `XSECTION.dat`** — el `inp.5` es idéntico |

### `Conjunto Simulaciones 2 - Análisis de cadenas`

CHAINS sobre los ocho isótopos del telurio natural, con ejecuciones monoisotópicas
para obtener las contribuciones y dos directorios auxiliares que generan las cintas
que el módulo lee.

### `Scripts`

Comprobaciones independientes, comentadas. **No usan la suite**: leen los mismos
ficheros de entrada que ACAB y resuelven el problema por otro camino.

---

## Nomenclatura de los casos

| Sufijo | Bibliotecas |
|---|---|
| `v.1 - Haffner` | Las del trabajo de referencia |
| `v.1b - Haffner con flujo térmico` | Ídem, con el convenio de flujo alternativo |
| `v.2 - DECAY (2007)` | Secciones de 211 grupos, decaimiento de 2007 |
| `v.3 - DECAY (2025)` | Secciones de 211 grupos, decaimiento de 2025 |
| `v.4 - DECAY-XSECTIONS (2025)` | Secciones de 175 grupos, decaimiento de 2025 |

Huellas SHA-256, primeros doce caracteres:

| Juego | `XSECTION.dat` | `DECAY.dat` | `REACTIONS.dat` |
|---|---|---|---|
| v.1 y v.1b | `20b528300f11` | `1eb135562e91` | `d8a55841cc56` |
| v.2 | `fc3667ee3428` | `17bbb43348e3` | `d8a55841cc56` |
| v.3 y referencia | `fc3667ee3428` | `8a610bac55ff` | `d8a55841cc56` |
| v.4 | `bc775556a87a` | `8a610bac55ff` | `03671ce7eb8b` |

**Los tres ficheros de datos nucleares son byte a byte los mismos en los cuatro
experimentos**: entre experimentos solo cambia el `inp.5`.

---

## Ficheros de cada carpeta

| Fichero | Qué es |
|---|---|
| `inp.5` | Entrada de ACAB |
| `fort.6` | Salida principal: inventarios y actividades |
| `XSECTION.dat`, `DECAY.dat`, `REACTIONS.dat` | Bibliotecas |
| `*.orig` | Versión anterior conservada, cuando la hubo |

**Nada se ha borrado.** Donde una simulación se reejecutó, la anterior se conserva
con extensión `.orig` y suele ser la evidencia de por qué se rehízo.

---

## Cómo reproducir un caso

```bash
cd "<carpeta del caso>"
acab.exe          # lee inp.5 y las tres bibliotecas del directorio
```

Y para comprobar el resultado sin usar la suite:

```bash
python Scripts/apoyo_ejemplo.py "<carpeta del caso>"
```

El acuerdo esperado con el `fort.6` es del orden del **0,04 %**, que es el suelo de
impresión del propio fichero de salida.

---

## Cinco cosas que hay que saber antes de reutilizar estos datos

Son las que costaron rondas de verificación durante el trabajo. **Léelas antes de
comparar cifras**, o las volverás a encontrar.

### 1. El `fort.6` imprime cuatro cifras significativas

Ninguna afirmación puede ser más fina que eso. Un cociente entre dos magnitudes
impresas tiene un suelo del orden de 10⁻³, y varias discrepancias que parecían reales
resultaron estar por debajo de él.

### 2. Las etiquetas de tiempo de las tablas van a tres cifras

Los encabezados de las tablas del `fort.6` redondean los instantes. **No los uses para
calcular**: los tiempos exactos están en el bloque 8 del `inp.5` y en su eco dentro
del propio `fort.6`. Dos ejecuciones con mallas distintas pueden imprimir la misma
etiqueta para instantes que difieren.

### 3. El origen temporal es el fin de la irradiación

Las tablas de enfriamiento cuentan desde el apagado, no desde el inicio de la
irradiación.

### 4. La malla de salida es también la malla de integración ⚠

**Es la advertencia importante.** Dos ejecuciones de ACAB con las tres bibliotecas
idénticas por huella, el mismo flujo y la misma composición, que difieran únicamente
en el número de instantes solicitados a la salida, **dan resultados distintos**: hasta
un **13,9 %** en el ¹³¹I.

El efecto aparece cuando λΔt del precursor rápido de la cadena —el ¹³¹Te, con 25 min—
cae en una banda acotada entre **6,8 y 31**, y se manifiesta como un exceso constante
por intervalo, de aproximadamente 0,96 veces la vida media del precursor. Afecta al
nucleido que **acumula**, no al que está en la banda.

**Las mallas de este repositorio están todas fuera de esa banda**, salvo las de
calibración, que existen precisamente para documentarlo. La referencia trabaja con
λΔt = 4,99, con un 33 % de margen.

*El mecanismo interno no está establecido*: dos hipótesis cayeron con contraejemplos.
El más limpio es que el ¹³¹Te (25,00 min) produce el efecto y el ¹²⁸I (24,98 min) no,
con el mismo λΔt.

### 5. La actividad impresa no es exactamente λN

El cociente entre la actividad que ACAB imprime y λ·N calculado con la semivida de su
propia biblioteca vale **1,0004 de media**, y es mayor que uno en las 23 especies
comprobadas. Es un desplazamiento sistemático en la conversión de átomos a actividad,
del orden del 0,04 %.

**Consecuencia práctica**: la tabla de átomos y la de actividades no dan exactamente
lo mismo al normalizar. Para cocientes entre casos, usa la tabla de átomos.

---

## Convenios

- **Actividad específica** referida a la masa de yodo, no a la del blanco, con la masa
  atómica de cada isótopo.
- **Techo sin portador**: λN_A/M para ¹³¹I puro, **4,6000×10⁹ MBq/g** con T½ = 693 200
  s y **4,5987×10⁹** con 693 400 s. Son dos y no deben unificarse.
- **Pureza radionucleídica** P = A(¹³¹I) / ΣA(isótopos de yodo), sin restricción de
  número másico.
- **Máximo de actividad**: se reporta por intervalo, no por instante, porque el suelo
  de impresión produce empates.

---

## Scripts

| Script | Qué hace |
|---|---|
| `apoyo_ejemplo.py` | Lee los datos nucleares de los ficheros de entrada, resuelve la cadena A = 131 en forma cerrada y contrasta con el `fort.6` |
| `control_bateman.py` | Control de linealidad de las ecuaciones de Bateman |

Están comentados y son cortos a propósito: **su valor está en ser independientes de
la suite y verificables a mano**, no en ser completos.

---

## Cita

Si utilizas este material, cita el TFG y el trabajo de referencia replicado:

> Haffner, R., Miller, W. H., & Morris, S. (2019). *Verification of I-131 Yield from
> the Neutron Irradiation of Tellurium.*

> Sanz, J., Cabellos, Ó., & García-Herranz, N. (2008). *ACAB — Activation Abacus
> Code.*

---

## Licencia y alcance

Los ficheros de entrada y salida de las simulaciones se publican para permitir la
verificación independiente de los resultados de la memoria. **Las bibliotecas de datos
nucleares están sujetas a las condiciones de sus evaluadores originales** y se
incluyen únicamente para hacer reproducible el cálculo.
