# MNOT

Its a library of solvers for non-obtuse triangulation problems.


## Getting started

### Installing dependencies

#### Unix / WSl
##### Ubuntu / Debian
```bash
sudo apt install -y build-essentail cmake
sudo apt install -y libcgal-dev
sudo apt install -y ipe
```
##### Arch
```bash
sudo pacman -S ipe cmake
sudo pacman -U file:///path/to/cgal-5.6.1-1-any.pkg.tar.zst.sig
```

Until CGAL 6.0 is supported it is necessary to download cgal-5.6.1-1 from the [Arch Linux Archive](https://archive.archlinux.org/packages/c/cgal/)



#### Windows without WSL
- Coming soon

### Building
```bash
git clone https://github.com/traMectory/MNOT.git
cd MNOT

mkdir build
cd build

cmake ../src
make
```

On Arch it might be neccesary to run `cmake -DCGAL_DIR=/usr/lib/cmake/CGAL ../src` instead.


## Usage
```bash
./mnot [path/to/instance_name.instance.json] [options]
```

The resulting triangulation will be saved in `.json` format according to [specification by CG:SHOP](https://cgshop.ibr.cs.tu-bs.de/competition/cg-shop-2025/#instance-format) in `solutions/instance_name.json`. If `-v` is specified, an `.ipe` file with the triangulation will be saved at `solutions/vis/instance_name.ipe`.

### Options
<table>
<tr>
<td><code>-v / --visualize</code></td><td>Show the resulting triangulation using IPE</td>
</tr>
<tr>
<td><code>-a ARG / --algorithm ARG</code></td><td>Use solver algorithm number <code>ARG</code></td>
</tr>
</table>

## Solvers
| Solver number  | Solver name |
| - | - |
| 0 | Stupid  |

## Adding your own solvers
- Coming soon