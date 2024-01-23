# Towards Adaptive Subspace Detection in Heterogeneous Environment

This GitHub repository was made to present the results of a two-step adaptive matched filter detector designed for a heterogeneous environment where the test data noise covariance may differ from the secondary set covariance in structure. In particular we use the model $\textbf{y}={\textbf{H}{\theta}}+{\textbf{B}\{\phi}}+\boldsymbol\xi$, where the hypothesis testing problem is a decision between ${\mathcal{H}_0}:
\textbf{y}=\textbf{B}{\phi}+{\xi} \sim \mathcal{C}\mathcal{N}(\textbf{B}{\phi},\sigma^2\textbf{R})$ with $\textbf{y}_k=\textbf{n}_k \sim \mathcal{C}\mathcal{N}(\textbf{0},\sigma^2_k\textbf{R}_s), k=1,\cdots,K$ and ${\mathcal{H}_1}:\textbf{y}=\textbf{H}{\theta}+\textbf{B}{\phi}+{\xi} \sim \mathcal{C}\mathcal{N}(\textbf{H}{\theta}+\textbf{B}{\phi},\sigma^2\textbf{R})$ with $\textbf{y}_k=\textbf{n}_k \sim \mathcal{C}\mathcal{N}(\textbf{0},\sigma^2_k\textbf{R}_s), k=1,\cdots,K$.

## Demo
+ Run the file "Main.m" for synthetic experiments.

## Results
The detection results (ROC) of the proposed method (red curve) vs AMF (with known covariance and estimated covariance) in different environments (HE,PHE, and HET):
![image](https://github.com/arekavandi/Heterogeneous_Detector/assets/101369948/7e996e74-f95b-42a6-9347-0c239a07547e)

## Citations
If you found this GitHub page helpful, please cite the following paper:

```
@article{rekavandi2024towards,
  title={Towards Adaptive Subspace Detection in Heterogeneous Environment},
  author={Miri Rekavandi, Aref},
  journal={arXiv},
  volume={},
  pages={1--8},
  year={2024},
  publisher={arXiv}
}
```
