# Towards Adaptive Subspace Detection in Heterogeneous Environment

This GitHub repository was made to present the results of an adaptive matched filter designed for a heterogeneous environment where the test data noise covariance may differ from the secondary set covariance in structure. In particular we use the model $\textbf{y}={\textbf{H}{\theta}}+{\textbf{B}\{\phi}}+\boldsymbol\xi$, where the hypothesis testing problem is a decision between ${\mathcal{H}_0}:
\textbf{y}=\textbf{B}{\phi}+{\xi} \sim \mathcal{C}\mathcal{N}(\textbf{B}{\phi},\sigma^2\textbf{R})$ with $\textbf{y}_k=\textbf{n}_k \sim \mathcal{C}\mathcal{N}(\textbf{0},\sigma^2_k\textbf{R}_s), k=1,\cdots,K$ and ${\mathcal{H}_1}:\textbf{y}=\textbf{H}{\theta}+\textbf{B}{\phi}+{\xi} \sim \mathcal{C}\mathcal{N}(\textbf{H}{\theta}+\textbf{B}{\phi},\sigma^2\textbf{R})$ with $\textbf{y}_k=\textbf{n}_k \sim \mathcal{C}\mathcal{N}(\textbf{0},\sigma^2_k\textbf{R}_s), k=1,\cdots,K$.

## Demo


## Results

## Citations
If you found this page helpful, please cite the following survey papers:

```
@article{rekavandi2021robust,
  title={Towards Adaptive Subspace Detection in Heterogeneous Environment},
  author={Miri Rekavandi, Aref},
  journal={arXiv},
  volume={30},
  pages={1--8},
  year={2024},
  publisher={arXiv}
}
```
