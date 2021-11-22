$$
\vec y_{n+1} = \vec y_n + h \sum_{s=1}^S b_s k_s
$$

$$
\begin{aligned}
\vec k_1 &= \vec 0 \\
\vec k_2 &= f \left( t_n + h c_1, ~ \vec y_n + h \sum_{i=1}^1 a_{1,i} \vec k_i \right) \\
\vec k_3 &= f \left( t_n + h c_2, ~ \vec y_n + h \sum_{i=1}^2 a_{2,i} \vec k_i \right) \\
\vec k_4 &= f \left( t_n + h c_3, ~ \vec y_n + h \sum_{i=1}^3 a_{3,i} \vec k_i \right) \\
\end{aligned}
$$

$$
\begin{aligned}
\delta t_s &= h c_s \\
\delta \vec y_s &= h \sum_{i=1}^{s} a_{s,i} \vec k_i
\end{aligned}
$$

$$
\vec k_{s+1} = f(t_n + \delta t_s, ~ \vec y_n + \delta \vec y_s )
$$
