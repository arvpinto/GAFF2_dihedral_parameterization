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

Then we extract the optimized geometries from the Gaussian scan output using the gaussian2xyz.py script (authored by Tomasz Borowski and Zuzanna Wojdyła) and run single-point calculations with a high level of theory:

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

#Finally we can use the gaussian_dihedral.py script to extract the energy profile for dihedral rotation:
python gaussian_dihedral.py 10 1 4 7 x*log > qm_scan.dat
</pre>

<div align="center">
    <img src="kernel_density_plot.png">
</div>

Note, while the dihedral angle scan is carried out with a lower level of theory, the single-point calculations should be carried out with an adequate method such as MP2/cc-pVTZ. Here we use M062X/6-311++G(d,p) to exemplify.

<br/>

Clean up the pc.pdb file to include only the PCA vectors:
<pre style="color: white; background-color: black;">
cat pc.pdb | head -n -2 | tail -n +6 | awk '{print $6,$7,$8}' > temp && mv temp pc.pdb
</pre>

<br/>
<h2> <p align="center"> <b>II - Clustering of PCA vectors and identification of representative frames</b> </p></h2>

<br/>

Now we run the <a href="https://arvpinto.github.io/3D_clustering_PCA/pca_dbscan_gmm.py" target="_blank">pca_dbscan_gmm.py</a> script to obtain the clusters and the representative frames.
The <a href="https://arvpinto.github.io/3D_clustering_PCA/pca_dbscan_gmm.py" target="_blank">pca_dbscan_gmm.py</a> script has the following usage:

<pre style="color: white; background-color: black;">
python pca_dbscan_gmm.py &lt;data_file&gt; &lt;eps&gt; &lt;min_samples&gt; &lt;n_components&gt;
</pre>

<p align="justify">The &lt;data_file&gt; should be the processed pc.pdb file, &lt;eps&gt; and &lt;min_samples&gt; define the parameters for outlier identification using the DBSCAN method, and &lt;n_components&gt; defines the number of clusters in the Gaussian Mixture Models clustering. The script produces a 3D plot of the PCA vectors, where the outliers are represented as black markers, the frames closest to the highest density points as white markers, and each cluster displays a different color. Additionally, the density distribution curves of each cluster are plotted against each PCA vector, with markers representing the identified frames.
Initially try different &lt;eps&gt; and &lt;min_samples&gt; values to see which and how many frames are being identified as outliers.
Once you have an adequate number of outliers, try different &lt;n_components&gt; values to identify which number of clusters is more suitable.
Also take a look at the kernel density plots to see if the density distributions have a regular shape, and the identified frames lie close to highest density points. </p>
<br/>

<pre style="color: white; background-color: black;">
Number of DBSCAN outliers: 29
Total number of clusters (GMM): 4
Cluster 0: 595 frames
Top 5 closest frames for Cluster 0: [ 578  721  681  647 1544]
Cluster 1: 1198 frames
Top 5 closest frames for Cluster 1: [1232 1380 1293 1919 1708]
Cluster 2: 463 frames
Top 5 closest frames for Cluster 2: [114  69  68  64  67]
Cluster 3: 215 frames
Top 5 closest frames for Cluster 3: [2015 2076 2050 2052 2054]
</pre>
<br>
<br>
<br>
<br>
<br>
<br>
<br>
<div style="display: flex; justify-content: center; align-items: center; height: 100vh;">
    <iframe src="https://arvpinto.github.io/3D_clustering_PCA/3d_plot.html" width="1904" height="894"></iframe>
</div>
<br>
<br>
<br>
<br>
<br>

<div align="center">
    <img src="kernel_density_plot.png">
</div>

A clusters.csv file is outputed with the cluster numbers that each frame corresponds to (outliers belong in the -1 cluster).
A frames.dat is ouputed with the top 5 frames that are closest to the highest density point of each cluster.

<br>
<h2> <p align="center"> <b>III - Frame extraction</b> </p></h2>

<br/>

Use the <a href="https://arvpinto.github.io/3D_clustering_PCA/extract_highdens.py" target="_blank">extract_highdens.py</a> script to extract the identified frames from the trajectory.
The <a href="https://arvpinto.github.io/3D_clustering_PCA/extract_highdens.py" target="_blank">extract_highdens.py</a> script usage follows:

<pre style="color: white; background-color: black;">
python extract_highdens.py &lt;xtc_file&gt; &lt;gro_file&gt; &lt;cluster_indices_file&gt; &lt;output_prefix&gt;
</pre>





