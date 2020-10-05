# Clustering-for-vectors-and-polygonal-curves
This project implements clustering algorithms for in Euclidian space R<sup>d</sup> and polygonal curves in Euclidian space R<sup>2</sup> using using the 8 combinations of the following variants. The Manhattan metric (L1) is used for vectors and the pseudo Dynamic Time Warping (DTW) for the curves.

**Initialization**
1. Random selection of K points/curves (simplest)
2. K-means++

**Assignment**
1. Lloyd's assignment (simplest approach)
2. Assignment by Range Search using LSH for vectors/curves (inverse assignment)

**Update**
1. Partitioning around medoids (PAM) \`a la Lloyds
2. Calculating Mean Vector/DTW Centroid Curve

## Installing

Download source code by typing:

```
https://github.com/PanPapag/Clustering-for-vectors-and-polygonal-curves.git
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details
