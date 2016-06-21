<h1> PageRank </h1>

aka <i>this was interesting.</i>

C implementation of the PageRank algorithm, with and without parallelization. Used as a group project for the <i>High Performance Computing</i> course held at Ca' Foscari University of Venice, master's degree in Computer science. The algorithm is implemented sequentially and then parallelized using the <i>openMP</i> library. 

<h1>The code</h1>

Several files are included:

<ul>
  <li><b>step1.c</b>, sequential implementation of PageRank. Uses a transposed adjacency matrix;</li>
  <li><b>step2.c</b>, sequential implementation of PageRank. Uses <i>compressed sparse row</i> organization of the adjacency matrix; </li>
  <li><b>step3.c</b>, parallel implementation of PageRank. Customize scheduling type (static, dynamic) and number of threads;</li>
  <li><b>step2mmap.c</b>, same as <i>step2.c</i> but using mmap;</li>
  <li><b>step3mmap.c</b>, same as <i>step3.c</i> but using mmap.</li>
</ul>

Compile the code using <a href="http://snap.stanford.edu/data/">these data sets</a> and see the magic. 

<h1>Acknowledgements</h1>

This project was done with with Gaia O. and Gianluca C..
