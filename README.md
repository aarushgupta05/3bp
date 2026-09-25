# Structure over Resolution: An Operating Envelope for Stable Integration of the Chaotic Three-Body Problem

**Author:** Aarush Gupta

This repository contains the code, numerical outputs, and figures used for the research paper **“Structure over Resolution: An Operating Envelope for Stable Integration of the Chaotic Three-Body Problem.”**

The study compares three numerical integrators for the three-body problem:

- Fixed-step leapfrog (LF)
- Time-symmetric adaptive leapfrog (TSALF)
- IAS15

They are tested across analytical three-body systems, binary-single systems, and the Sun-Earth-Moon system. The aim is to identify which integrator is most suitable for different dynamical regimes, accuracy requirements, and computational costs.

## Repository Structure

- `src/` — core N-body simulation and integrator code
- `scripts/` — experiment, analysis, and figure-generation scripts
- `outputs/` — numerical results and generated data tables
- `figures/` — figures used in the paper
- `REPRODUCE.md` — detailed instructions for reproducing the results
- `requirements.txt` — Python dependencies
- `CITATION.cff` — citation information
- `LICENSE` — MIT License

## Installation

The simulations were run using Python 3.11.

Create a Python environment and install the required packages:

```bash
pip install -r requirements.txt
```

## Reproducing the Results

Detailed reproduction instructions are provided in [`REPRODUCE.md`](REPRODUCE.md).

The paper figures can be regenerated using:

```bash
python scripts/make_fig1.py
python scripts/make_fig2.py
python scripts/make_fig3.py
python scripts/make_fig4.py
```

## Main Result

The results show that there is no single best integrator for every three-body system. Fixed-step leapfrog can be efficient when uniform resolution is sufficient, TSALF is useful when the required resolution changes during the trajectory, and IAS15 performs well for tightly bound systems requiring high accuracy.

The Sun-Earth-Moon experiment also shows that good global energy conservation does not necessarily guarantee that the internal physical structure of a simulation is correct, so additional diagnostics such as pair separation are important.

## Citation

If you use this code or results, please use the citation information provided in [`CITATION.cff`](CITATION.cff).

## License

This project is released under the MIT License. See [`LICENSE`](LICENSE) for details.
