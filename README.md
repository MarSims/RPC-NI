# RPC-NI

Representative Point-Based Clustering with Neighborhood Information for Complex Data Structures (RPC-NI)

*Zhongju Shang, Yaoguo Dang, Haowei Wang, Sifeng Liu*

Published in: IEEE Transactions on Cybernetics ( Volume: 55, Issue: 4, April 2025)

https://ieeexplore.ieee.org/document/10887271

![Image text](https://github.com/MarSims/img-folder/blob/main/framework-R1.png)

# Abstract

Discovering clusters remains challenging when dealing with complex data structures, including those with varying densities, arbitrary shapes, weak separability, or the presence of noise. In this article, we propose a novel clustering algorithm called representative point-based clustering with neighborhood information (RPC-NI), which highlights the significance of neighborhood information often overlooked by existing clustering methods. The proposed algorithm first introduces a new local centrality metric that integrates both neighborhood density and topological convergence to identify core representative points, effectively capturing the structural characteristics of the data. Subsequently, a density-adaptive distance is defined to evaluate dissimilarities between these core representative points, and such distance is used to construct a minimum spanning tree (MST) over these points. Finally, an MST-based clustering algorithm is employed to yield the desired clusters. Incorporating neighborhood information enables RPC-NI to comprehensively determine representative points, and having multiple representative points per cluster allows RPC-NI to adapt to clusters of arbitrary shapes, varying densities, and different sizes. Extensive experiments on widely used datasets demonstrate that RPC-NI outperforms baseline algorithms in terms of clustering accuracy and robustness. These results provide further evidence for the importance of incorporating neighborhood information discovering clusters with complex structures.

# 

```
@article{shang2025representative,
  title={Representative Point-Based Clustering With Neighborhood Information for Complex Data Structures},
  author={Shang, Zhongju and Dang, Yaoguo and Wang, Haowei and Liu, Sifeng},
  journal={IEEE Transactions on Cybernetics},
  volume={55},
  number={4},
  pages={1620 - 1633},
  year={2025},
  publisher={IEEE} 
}
```

# 

This code has been evaluated on Matlab 2023b

Email: shangzhongju@163.com
