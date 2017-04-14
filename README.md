# Monte-Carlo study of the 2D Ising model
This project deals with the thermodynamics of the ferromagnetic Ising Model - a 2D array of spins with simple nearest neighbour interactions. This model is of interest because it displays a phase transition from disorder → order at a critical temperature of β<sub>c</sub>=0.4407. This comes from the exact solution in the large lattice limit by L. Onsager.

## Build and run
Make sure you have cpgplot installed, the C-callable version of pgplot <http://www.astro.caltech.edu/~tjp/pgplot/>

### Command line
Build:
```
make proj_4
```

Then run:
```
./proj_4
```

### Within CLion, JetBrains' IDE for C projects
The `CMakeLists.txt` file should give you everything you need to specify build targets and required libraries, allowing you to build and run within the IDE.

## Preparation of report
You can build the finished pdf of the report using the `Project4.tex` file and your favourite LaTeX compilation program. I used TexStudio on Mac <http://www.texstudio.org>