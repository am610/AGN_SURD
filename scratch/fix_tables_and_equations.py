import re

# Load V3.tex
with open('/Users/ayan/Programs/SURD/overleaf_draft/V3.tex', 'r') as f:
    text = f.read()

# 1. Replace figure* float placement arguments completely
# (e.g. \begin{figure*}[t] -> \begin{figure*})
# Let's use regex for safety:
text = re.sub(r'\\begin\{figure\*\}\[[a-zA-Z!]+\]', r'\\begin{figure*}', text)

# 2. Replace all \begin{table}[ht] with \begin{table*} and \end{table} with \end{table*}
# We can do this with simple replace calls:
text = text.replace(r'\begin{table}[ht]', r'\begin{table*}')
text = text.replace(r'\end{table}', r'\end{table*}')

# 3. Replace Equation 24 with split version
old_eq24 = r"""\begin{equation}
\mathbf{X}(t)=
\Bigl[
Q_1(t),Q_1(t-\Delta T_1),\ldots,Q_1(t-\Delta T_p),\;
Q_2(t),Q_2(t-\Delta T_1),\ldots,Q_N(t-\Delta T_p)
\Bigr].
\label{eq:lagged_vector_general}
\end{equation}"""

new_eq24 = r"""\begin{equation}
\begin{split}
\mathbf{X}(t) = \Bigl[ & Q_1(t), Q_1(t-\Delta T_1), \ldots, Q_1(t-\Delta T_p), \\
                       & Q_2(t), Q_2(t-\Delta T_1), \ldots, Q_N(t-\Delta T_p) \Bigr].
\end{split}
\label{eq:lagged_vector_general}
\end{equation}"""

text = text.replace(old_eq24, new_eq24)

# 4. Replace Equation 26 with split version
old_eq26 = r"""\begin{equation}
\mathbf{X}(t)=
\bigl[
\text{past lags},
\text{present variables},
\text{contemporaneous non-target variables}
\bigr].
\end{equation}"""

new_eq26 = r"""\begin{equation}
\begin{split}
\mathbf{X}(t) = \bigl[ & \text{past lags}, \text{ present variables}, \\
                       & \text{contemporaneous non-target variables} \bigr].
\end{split}
\end{equation}"""

text = text.replace(old_eq26, new_eq26)

# Overwrite V3.tex
with open('/Users/ayan/Programs/SURD/overleaf_draft/V3.tex', 'w') as f:
    f.write(text)

print("Finished adjusting tables, figures, and equations in V3.tex")
