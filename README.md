<p align="justify"><b>This tutorial explains how to reparameterize molecular mechanics (MM) dihedral angles using a quantum mechanical (QM) approach.</b></p>

<p align="justify">
It requires Gaussian, GROMACS and a recent version of AmberTools.
</p>

---

<br>
<h2> <p align="center"> <b>I - Dihedral angle scan with QM </b> </p></h2>

<br/>

We start by creating our molecule in GaussView and performing a dihedral scan using the ModRedundant input section: 
<pre style="color: white; background-color: black;">
# opt=modredundant b3lyp/6-31g(d,p) scf=tight

Title Card Required

0 1
 C                  1.39817300    0.52720100    0.22441400
 H                  2.01440700    1.35294300   -0.16438900
 H                  1.44317500    0.57662300    1.32520700
 C                 -0.04224600    0.70264600   -0.23892700
 H                 -0.06319300    0.64168000   -1.33558100
 H                 -0.39046500    1.70533500    0.03570800
 C                 -1.00011300   -0.33618100    0.34929500
 H                 -0.65273900   -1.34722800    0.09797700
 H                 -1.00776200   -0.25692900    1.44251900
 O                  1.86718200   -0.73237700   -0.25133700
 H                  2.76112600   -0.86758100    0.08668800
 O                 -2.34813900   -0.12614100   -0.05961100
 H                 -2.39178000   -0.29869800   -1.00922800

10 1 4 7 S 60 6.0
12 7 4 1 F

</pre>

Here we use the 1,3-propanediol molecule as an example, scanning the 10 1 4 7 dihedral while keeping the 12 7 4 1 dihedral frozen to avoid energy jumps. Note that the initial geometry of the scan should correspond to a global minimum in the potential energy surface. If possible, use a higher level of theory.

<br/>

Then we extract the optimized geometries from the Gaussian scan output using the gaussian2xyz.py <a href="https://arvpinto.github.io/3D_clustering_PCA/pca_dbscan_gmm.py" target="_blank">pca_dbscan_gmm.py</a> script (authored by Tomasz Borowski and Zuzanna Wojdyła) and run single-point calculations with a high level of theory:

<pre style="color: white; background-color: black;">
#The script generates the scan_geoms.xyz file which has all the coordinates and corresponding energies
python gaussian2xyz.py propanediol_scan.log scan > scan_geoms.xyz

#We split the file to create individual coordinate files
split -l 15 scan_geoms.xyz

#Using a for loop, we create new Gaussian input files, which will be used for single-point calculations
for i in x*; do 
   (echo "#p M062X/6-311++G(d,p)" ; echo "" ; echo "SP" ; echo "" ; echo "0 1" ; tail -n +3 "$i" ; echo "") > "$i".com 
done

#And we run the calculations in the background:
nohup $(for i in *com ; do g09 "$i" ; done) &

#Finally we can use the gaussian_dihedral.py <a href="https://arvpinto.github.io/3D_clustering_PCA/pca_dbscan_gmm.py" target="_blank">pca_dbscan_gmm.py</a> script to extract the energy profile for dihedral rotation:
python gaussian_dihedral.py 10 1 4 7 x*log > qm_scan.dat
</pre>

<div align="center">
    <img src="kernel_density_plot.png">
</div>

Note, while the dihedral angle scan is carried out with a lower level of theory, the single-point calculations should be carried out with an adequate method such as MP2/cc-pVTZ. Here we use M062X/6-311++G(d,p) to exemplify.

<br/>
<h2> <p align="center"> <b>II - Dihedral angle scan with MM </b> </p></h2>

<br/>

To perform an MM dihedral scan, the molecule is first parameterized with GAFF2:
<pre style="color: white; background-color: black;">
antechamber -fi gout -i propanediol.log -fo mol2 -o propanediol.mol2 -nc 0 -c abcg2 -pf y -at gaff2
antechamber -fi mol2 -i propanediol.mol2 -fo ac -o propanediol.ac -pf y
prepgen -i propanediol.ac -o propanediol.prepin
parmchk2 -i propanediol.prepin -f mol2 -o propanediol.frcmod -s gaff2

tleap
>source leaprc.gaff2
>loadamberprep propanediol.prepin
>mol = loadmol2 propanediol.mol2
>saveoff mol propanediol.lib
>quit

tleap
>source leaprc.gaff2
>loadamberparams propanediol.frcmod
>loadoff propanediol.lib
>loadamberprep propanediol.prepin
>mol = loadmol2 propanediol.mol2
>saveAmberParm mol propanediol.prmtop propanediol.rst7
>savepdb mol propanediol.pdb
>quit
</pre>

Then the coordinates and parameters are converted to GROMACS format:
<pre style="color: white; background-color: black;">
python amb2gmx.py propanediol.prmtop propanediol.rst7
gmx editconf -f propanediol_converted.gro -bt triclinic -d 1.0 -o propanediol_converted.gro
</pre>

To carry out the dihedral scan, we introduce the following section in the topology after the parameters:
<pre style="color: white; background-color: black;">
#ifdef POSRES
[ dihedral_restraints ]
10 1 4 7 1 DIHE_VALUE 0 10000
12 7 4 1 1 177.95895 0 10000
#endif
</pre>

And then we run the calculations with the min_steep_restr.mdp <a href="https://arvpinto.github.io/3D_clustering_PCA/pca_dbscan_gmm.py" target="_blank">pca_dbscan_gmm.py</a> file:
<pre style="color: white; background-color: black;">
for i in $(tail -n +2 qm_scan.dat | awk '{print $1}'); do 
    cp propanediol_converted.top propanediol_converted_dihe.top 
    sed -i 's/DIHE_VALUE/'"$i"'/g' propanediol_converted_dihe.top 
    gmx grompp -f min_steep.mdp -c propanediol_converted.gro -p propanediol_converted_dihe.top -o dihe_"$i".tpr 
    gmx mdrun -deffnm dihe_"$i"  
done
</pre>

This will result in a dihedral energy profile with the dihedral term that we aim to parameterize included, however we need to calculate the energy profile without this dihedral term (see https://pubs.acs.org/doi/10.1021/acs.jpca.0c10845). This can be done by deleting the term from the topology file and running single-point calculations on the previously produced structures with the min_steep_sp.mdp <a href="https://arvpinto.github.io/3D_clustering_PCA/pca_dbscan_gmm.py" target="_blank">pca_dbscan_gmm.py</a> file:
<pre style="color: white; background-color: black;">
for i in dihe_*gro; do 
   gmx grompp -f min_steep_sp.mdp -c "$i" -p propanediol_converted.top -o "$(echo "$i" | sed 's/\.gro//')".tpr 
   gmx mdrun -deffnm "$(echo "$i" | sed 's/\.gro//')" 
done
</pre>

Then the zero-torsion energy profile can be extracted with the gromacs_dihedral.py <a href="https://arvpinto.github.io/3D_clustering_PCA/pca_dbscan_gmm.py" target="_blank">pca_dbscan_gmm.py</a> script:
<pre style="color: white; background-color: black;">
python3 gromacs_dihedral.py dihe_*log
</pre>

<div align="center">
    <img src="kernel_density_plot.png">
</div>

<br/>
<h2> <p align="center"> <b>III - Fitting the dihedral energy term </b> </p></h2>

Now we can use the least_squares_fit.py script to fit Fourier (cos/sin) series to the target dihedral energy profile:
<pre style="color: white; background-color: black;">
python fit_dihedral_gromacs.py --qm qm_scan.dat --mm mm_scan.dat --out dihedral.itp --nmax 1 --refine
</pre>
Here we only have one term associated with the dihedral, so we use --nmax 1.

Finally, we can replace the term in the original topology and re-run the MM single-point calculations with this term included to see if the fit adequately leads to the reproduction of the QM torsional profile:

