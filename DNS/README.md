# Drop impact onto a moving pool

Direct numerical simulation code infrastructure for **three-dimensional drop impact onto a moving pool**, supporting collaborative work with the Oxford Fluids Lab at Oxford and the Harris Lab at Brown. The code complements the [associated arXiv preprint](https://arxiv.org/abs/2511.03682), and provides tools for parameter specification, verification and visualisation.

<img width="100%" alt="3D_SimulationSnapshot" src="https://github.com/user-attachments/assets/4b56ec1d-5d7d-4acc-b828-ccb7533f9303" />

---

## 📌 Features

✅ Three-dimensional Navier-Stokes solver for drop-pool impact scenarios   
✅ Two-phase, coalescing fluid volume implementation suitable for splashing dynamics
✅ High-resolution output and animation capabilities  
✅ Supplementary movies for publication and presentation  

---

## 🛠️ Installation

### 1. Requirements

- [Basilisk](http://basilisk.fr/) (compiled with `qcc`)
- C compiler
- Visualization tools (e.g., ffmpeg, imagemagick, gnuplot)
- Recommended:  
  ```bash
  sudo apt install ffmpeg imagemagick gnuplot
  ```

### 2. Clone the repository

```bash
git clone https://github.com/rcsc-group/MovingPoolImpact
cd MovingPoolImpact/DNS
```

### 3. Install Basilisk

See Basilisk’s [installation page](http://basilisk.fr/src/INSTALL) for instructions.

### 4. Running the code

After Basilisk is set up, run the driver code using the provided shell script:

```bash
sh run_movingpool.sh
```

The code organizes outputs into folders with summaries and movies.

---

## ⚙️ Key Simulation Options

Edit parameters across various layers to control:

- Maximum resolution level
- Impingement angle
- Initial drop velocity
- Horizontal pool velocity
- Drop radius
- Pool Depth
- Computational Domain size
- Simulation end time
- Output frequency and visualization options (inside the driver code itself)

---

## 📁 Folder Structure

Kept simple, just a driver code, a running shell script and the current README.

---

## 📊 Outputs

Generates:

Summary files and movies in organized output folders:
- `.mp4` animations for simulation data
- Data files (interface data coordinates) for further analysis

Visualisation can be toggled depending on your architecture and needs.

---

## 📚 Citation

If you use this code or data in your work, please cite the associated preprint (for now):

> Thomas C. Sykes, Luke F.L. Alventosa, J. Rafael Castrejon-Pita, Radu Cimpeanu, Daniel M. Harris, Alfonso A. Castrejon-Pita (2025). Fast droplet impact onto slowly moving deep pools. arXiv:2511.03682.

BibTeX:
```bibtex
@article{Sykes2025arxiv,
 author = {Sykes, Thomas C. and Alventosa, Luke F.L. and Castrej{\'o}n-Pita, J. Rafael and Cimpeanu, Radu and Harris, Daniel M. and Castrej{\'o}n-Pita, Alfonso A.},
 title = {Fast droplet impact onto slowly moving deep pools},
 journal = {arXiv preprint arXiv:2511.03682},
 doi = {https://doi.org/10.48550/arXiv.2511.03682},
 volume = {},
 pages = {},
 year = {2025}
}
```

---

## 🧑 Contributing

Feel free to:

- Fork this repo
- Open issues
- Submit pull requests

---

```
