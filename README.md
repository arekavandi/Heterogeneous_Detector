[![arXiv](https://img.shields.io/badge/arXiv-2306.04670-b31b1b.svg)](https://arxiv.org/abs/2401.12469)
# Towards Adaptive Subspace Detection in Heterogeneous Environments

This GitHub repository was made to present the results of a two-step adaptive matched filter detector designed for a heterogeneous environment where the test data noise covariance may differ from the secondary set covariance in structure. In particular we use the model $\textbf{y}={\textbf{H}{\theta}}+{\textbf{B}\{\phi}}+\boldsymbol\xi$, where the hypothesis testing problem is a decision between ${\mathcal{H}_0}:
\textbf{y}=\textbf{B}{\phi}+{\xi} \sim \mathcal{C}\mathcal{N}(\textbf{B}{\phi},\sigma^2\textbf{R})$ with $\textbf{y}_k=\textbf{n}_k \sim \mathcal{C}\mathcal{N}(\textbf{0},\sigma^2_k\textbf{R}_s), k=1,\cdots,K$ and ${\mathcal{H}_1}:\textbf{y}=\textbf{H}{\theta}+\textbf{B}{\phi}+{\xi} \sim \mathcal{C}\mathcal{N}(\textbf{H}{\theta}+\textbf{B}{\phi},\sigma^2\textbf{R})$ with $\textbf{y}_k=\textbf{n}_k \sim \mathcal{C}\mathcal{N}(\textbf{0},\sigma^2_k\textbf{R}_s), k=1,\cdots,K$.

## Demo
+ Run the file "Main_new.m" for synthetic experiments.

## Results
The detection results (ROC) of the proposed method (black curve) vs AMF (with known covariance and estimated covariance) and ASD (with known covariance and estimated covariance) in different environments (HE,PHE, and HET):
![image](https://github.com/arekavandi/Heterogeneous_Detector/assets/101369948/65b37abd-a5e7-44ab-adb3-04a37a763a07)



## Citations
If you found this GitHub page helpful, please cite the following paper:

```
@article{rekavandi2024towards,
  title={Towards Adaptive Subspace Detection in Heterogeneous Environments},
  author={Rekavandi, Aref Miri},
  journal={arXiv preprint arXiv:2401.12469},
  year={2024}
}
```
