$
z_{k+1}=z_k-\frac{\eta{i_k}\Delta{t}}{Q}+w_k\\
v_k=OCV{(z_k)}-i_kR+v_k\\$<br>

$
\begin{array}{lccccc}
\hline\rule{0pt}{3.5ex} \text { Method } & \gamma & \alpha_0^{(\mathrm{m})} & \alpha_k^{(\mathrm{m})} & \alpha_0^{(\mathrm{c})} & \alpha_k^{(\mathrm{c})} \\
\text { CDKF } & h & \frac{h^2-L}{h^2} & \frac{1}{2 h^2} & \frac{h^2-L}{h^2} & \frac{1}{2 h^2}\rule[-2ex]{0pt}{5ex}\\
\hline
\end{array}
\\$
$_{\quad {}^*h \text{ may take any positive value. For Gaussian RVs, } h = \sqrt{3}}\\\quad{_{L\text{ is the 
dimension of state space vector }x}}$
### <b>SPKF step 1a:</b> State estimate time update. 

- Augmented aposteriori state estimate vector for previous time interval :

$\quad \hat{\chi}{^{a,+}_{k-1}}=\begin{bmatrix}{\hat{\chi}{^+_{k-1}}},&{\bar{w}},&{\bar{v}}\end{bmatrix}$<br><br>
- Augmented aposteriori state covariance estimate vector for previous time interval:

$\quad\sum{^{a,+}_{\tilde{x},k-1}}=\mathrm{diag}({\sum{^{+}_{\tilde{x},k-1}}},\sum{\tilde{w}},\sum{\tilde{v}})$<br><br>
- To generate the p+1 augmented sigma points

$\quad\chi{^{a,+}_{k-1}}=\left\{{\hat{\chi}{^{a,+}_{k-1}}},\quad\hat{\chi}{^{a,+}_{k-1}}+\gamma\sqrt{\sum
{^{a,+}_{\tilde{x},k-1}}},\quad\hat{\chi}{^{a,+}_{k-1}}-\gamma\sqrt{\sum{^{a,+}_{\tilde{x},k-1}}}\right\}\\ \quad \rule[-1.5ex]{0pt}{4.5ex}
\mathcal{X}_{k, i}^{x,-}=f\left(\mathcal{X}_{k-1, i}^{x,+}, u_{k-1}, \mathcal{X}_{k-1, i}^{w,+}\right)
$
$
\begin{aligned}
\quad\hat{x}_k^{-} & =\mathbb{E}\left[f\left(x_{k-1}, u_{k-1}, w_{k-1}\right) \mid \mathbb{Y}_{k-1}\right]
 \approx \sum_{i=0}^p \alpha_i^{(\mathrm{m})} f\left(\mathcal{X}_{k-1, i}^{x,+}, u_{k-1}, \mathcal{X}_
 {k-1, i}^{w,+}\right)
 =\sum_{i=0}^p \alpha_i^{(\mathrm{m})} \mathcal{X}_{k, i}^{x,-}
\end{aligned}
$
### <b>SPKF step 1b:</b> Error covariance time update.

- Using the *apriori* sigma points from step 1a, the *apriori* covariance
estimate is computed as :

$\quad\large{{\Sigma}^{-}_{\tilde{x},k} = \sum_{i=0}^{p} \alpha_i^{(c)} (\mathcal{X}_{k,i}^{x,-} - \hat{x}
_k^{-}) (\mathcal{X}_{k,i}^{x,-} - \hat{x}_k^{-})^T
}$
### <b>SPKF step 1c:</b> Estimate system output $y_k$.
- First, we compute the points :

$\quad
\mathcal{Y}_{k,i} = h(\mathcal{X}_{k,i}^{x,-}, u_k, \mathcal{X}_{k-1,i}^{v,+})
$
- The output estimate is then

$\quad
\hat{y}_k = \mathbb{E}[h(x_k, u_k, v_k) \mid \mathbb{Y}_{k-1}] \approx \sum_{i=0}^{p} \alpha_i^{(m)} h(\mathcal{X}_{k,i}^x, u_k, \mathcal{X}_{k-1,i}^{v,+}) = \sum_{i=0}^{p} \alpha_i^{(m)} \mathcal{Y}_{k,i}.
$
### <b>SPKF step 2a:</b> Estimator gain matrix $L_k$
- To compute the estimator gain matrix, we must first compute the
required covariance matrices.

$\large{\quad
\Sigma_{\tilde{y} , k}=\sum_{i=0}^{p} \alpha_{i}^{( \mathrm{c} )} \big( \mathcal{Y}_{k , i}-\hat{y}_{k} \big) \big( \mathcal{Y}_{k , i}-\hat{y}_{k} \big)^{T}}\\
\rule[-.5ex]{0pt}{3ex}
\large{\quad
\Sigma_{\tilde{x} \tilde{y} , k}^{-}=\sum_{i=0}^{p} \alpha_{i}^{( \mathrm{c} )} \big( \mathcal{X}_{k , i}^{x ,-}-\hat{x}_{k}^{-} \big) \big( \mathcal{Y}_{k , i}-\hat{y}_{k} \big)^{T}}$

- Then, we simply compute $L_k=\sum^-_{\tilde{x}\tilde{y},k}\sum^{-1}_{\tilde{y},k}$

### <b>SPKF step 2b:</b> State estimate measurement update.
- The state estimate is computed as

$\quad\hat{x}^+_k=\hat{x}^-_k+L_k{(y_k-\hat{y}_k)}$

### <b>SPKF step 2c:</b> Error covariance measurement update.
- The final step is calculated directly from the optimal formulation:

$\large\quad\sum{^{+}_{\tilde{x},k}}=\sum_{\tilde{x},k}^--
L_k\sum_{\tilde{y},k}L_k^T$
