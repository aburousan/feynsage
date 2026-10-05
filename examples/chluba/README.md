# QED results of Chluba's thesis, recomputed

J. Chluba, *Spectral distortions of the cosmic microwave background* (2005) uses Compton and double
Compton scattering and meets the infrared divergence of double Compton emission. Here each QED result
is computed from the Feynman diagrams twice: with feynsage + FORM (Jupyter notebooks) and with FeynCalc
(`*.wl`, FeynHelpers/Package-X for the loop).

- `compton_double_compton.ipynb`: Compton scattering (Peskin (5.87), Klein-Nishina, Thomson), double
  Compton with its six diagrams (Mandl and Skyrme's X) and the soft limit
- `ir_cancellation.ipynb`: the infrared pole of the vertex against the soft-photon pole, F2 = alpha/2pi,
  Lightman's law
- `double_compton.py`: the bookkeeping of the six double Compton diagrams, used by the first notebook
  and by section 18 of `../feynsage_walkthrough.ipynb`

Open the notebooks with the SageMath kernel. The FeynCalc versions run from the repository root:

    wolframscript -file examples/chluba/compton_dc_feyncalc.wl
    wolframscript -file examples/chluba/ir_cancellation_feyncalc.wl

The `.wl` files need FeynCalc and FeynHelpers in the Wolfram kernel that runs them. If `wolframscript`
picks another kernel (the output then shows unevaluated `DiracTrace[...]`), run them with the kernel
of the Mathematica that has FeynCalc, for example on macOS
`/Applications/Wolfram.app/Contents/MacOS/WolframKernel -noprompt -script examples/chluba/compton_dc_feyncalc.wl`.

| Result | feynsage + FORM | FeynCalc |
|---|---|---|
| Compton: sum of \|M\|^2 = 8 e^4 [...] (textbook) | exact | exact |
| Klein-Nishina and sigma_T = 8 pi r0^2/3 | exact | exact |
| double Compton, 6 diagrams: sum of \|M\|^2 = 4 e^6 X (Mandl & Skyrme, thesis eq. D.1) | exact at 4 points | exact at 4 points |
| soft limit: sum \|M_DC\|^2 -> e^2 S(K2) sum \|M_C\|^2 | ratio 1 - 3e-7 at K2 x 1e-8 | |
| Lightman's 4 alpha/(3 pi) (thesis eq. 4.24) from the eikonal factor | exact | |
| IR pole of the vertex = soft-photon emission pole (dimensional regularisation) | 1e-16 | 1e-15 |
| F2 -> alpha/(2 pi) | yes | yes |

So the thesis's |M|^2 = e^6 X is the average over the 2 x 2 initial states, and its lowest-frequency
cutoff is what dimensional regularisation (or a photon mass) turns into a log of the energy resolution:
the infrared divergence cancels against the one-loop virtual correction (Bloch-Nordsieck).
The FeynCalc outputs go to `examples/chluba/out/`.
