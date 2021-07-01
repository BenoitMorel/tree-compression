# Tree compression

Simple prototype implementation of tree compression methods (a simple compression and a compression using the RF metric) in C / C++.

## Getting Started

(In this fork, one can also run ./install.sh to download and build the dependencies, and build the project)


Download the repository, install the required libraries PLL modules and SLSL-lite. 
Once you have installed the two libraries, you can compile the project running
```
make
```
You then can run a simple compression on two trees with 500 taxa stored in the folder 500/ by running
```
./main ../data/500/tree_1.nwk ../data/500/tree_2.nwk
```
Also, you can change the trees by changing the number of the tree (e.g. using tree_10.nwk). There are 200 trees in the folder. Additionally, trees with 27 taxa are stored in 027/. 

### Prerequisites

To be able to run the tree compression, you will need to download and install the PLL modules 
```
git clone https://github.com/ddarriba/pll-modules.git
```
as well as the SLSL library
```
git clone https://github.com/simongog/sdsl-lite.git
```

## Authors

* **Axel Trefzer** - [axeltref](https://github.com/axeltref)
* **Alexis Stamatakis** - *Ideas and advice* - [stamatak](https://github.com/stamatak)
* **Simon Gog** - *Initial ideas* - [simongog](https://github.com/simongog)

