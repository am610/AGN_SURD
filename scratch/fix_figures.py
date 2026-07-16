# Load V3.tex
with open('/Users/ayan/Programs/SURD/overleaf_draft/V3.tex', 'r') as f:
    text = f.read()

# Define a replacement mapping for the figures
# We want to change \begin{figure}[ht] (or \begin{figure}[h]) to \begin{figure*}[t]
# and the matching \end{figure} to \end{figure*}

# Let's inspect the files around each figure and do a precise regex swap.
# In LaTeX, figure blocks look like:
# \begin{figure}[ht]
# ...
# \end{figure}
# We can find all of them and replace them.

import re

# We search for \begin{figure}[...] or \begin{figure} and replace with \begin{figure*}[t]
# and the matching \end{figure} with \end{figure*}.
# Since we want to make sure we don't accidentally mess up formatting, let's do a find-and-replace
# of the specific blocks:

old_fig1 = r"""\begin{figure}[ht]
\centering
\includegraphics[width=\textwidth]{figure1_light_curves.png}
\caption{Top panel: Standardized 5100~\AA\ continuum light curve for NGC~5548 over the strictly aligned overlap window (MJD 47512--49255). Bottom panel: Standardized velocity-resolved H\(\beta\) components (Blue Wing, Line Core, and Red Wing) resampled onto the shared 1.0-day grid.}
\label{fig:light_curves}
\end{figure}"""

new_fig1 = r"""\begin{figure*}[t]
\centering
\includegraphics[width=\textwidth]{figure1_light_curves.png}
\caption{Top panel: Standardized 5100~\AA\ continuum light curve for NGC~5548 over the strictly aligned overlap window (MJD 47512--49255). Bottom panel: Standardized velocity-resolved H\(\beta\) components (Blue Wing, Line Core, and Red Wing) resampled onto the shared 1.0-day grid.}
\label{fig:light_curves}
\end{figure*}"""

old_fig2 = r"""\begin{figure}[ht]
\centering
\includegraphics[width=\textwidth]{figure2_surd_lag_scans.png}
\caption{Decomposition of normalized predictive information (left column, where components represent the fraction of the joint mutual information $I(Y; X_1, X_2)$ and sum to $1.0$ at every lag: Synergy $\widehat{S}_{12}$ in solid purple, Redundancy $\widehat{R}_{12}$ in dashed gray, Unique Continuum $\widehat{U}_1$ in dotted blue, and Unique Wing/Core $\widehat{U}_2$ in dash-dotted red) and normalized information leak $\mathcal{L} = H(Y \mid X_1, X_2) / H(Y)$ (right column) as a function of lag up to 120 days for the three H\(\beta\) target components: Core (top row), Red Wing (middle row), and Blue Wing (bottom row). Solid curves represent the Monte Carlo median, and the dark and light shaded envelopes denote the $1\sigma$ (16th--84th percentiles) and $2\sigma$ (2.5th--97.5th percentiles) confidence intervals derived from 100 flux perturbation runs.}
\label{fig:surd_lag_scans}
\end{figure}"""

new_fig2 = r"""\begin{figure*}[t]
\centering
\includegraphics[width=0.95\textwidth]{figure2_surd_lag_scans.png}
\caption{Decomposition of normalized predictive information (left column, where components represent the fraction of the joint mutual information $I(Y; X_1, X_2)$ and sum to $1.0$ at every lag: Synergy $\widehat{S}_{12}$ in solid purple, Redundancy $\widehat{R}_{12}$ in dashed gray, Unique Continuum $\widehat{U}_1$ in dotted blue, and Unique Wing/Core $\widehat{U}_2$ in dash-dotted red) and normalized information leak $\mathcal{L} = H(Y \mid X_1, X_2) / H(Y)$ (right column) as a function of lag up to 120 days for the three H\(\beta\) target components: Core (top row), Red Wing (middle row), and Blue Wing (bottom row). Solid curves represent the Monte Carlo median, and the dark and light shaded envelopes denote the $1\sigma$ (16th--84th percentiles) and $2\sigma$ (2.5th--97.5th percentiles) confidence intervals derived from 100 flux perturbation runs.}
\label{fig:surd_lag_scans}
\end{figure*}"""

old_fig3 = r"""\begin{figure}[ht]
\centering
\includegraphics[width=\textwidth]{figure3_robustness_and_nulls.png}
\caption{Panel A: Sensitivity of Core H\(\beta\) synergy to the number of histogram bins ($n_{\mathrm{bins}} = 4, 6, 8, 10, 12$) over lags 1--60 days, showing strong qualitative stability. Panel B: Real synergy curve for Core H\(\beta\) compared with the median and 95\% percentile envelopes for circular-shift (blue) and block-bootstrap (red) surrogate runs, indicating local significance relative to standard nulls.}
\label{fig:robustness_and_nulls}
\end{figure}"""

new_fig3 = r"""\begin{figure*}[t]
\centering
\includegraphics[width=0.95\textwidth]{figure3_robustness_and_nulls.png}
\caption{Panel A: Sensitivity of Core H\(\beta\) synergy to the number of histogram bins ($n_{\mathrm{bins}} = 4, 6, 8, 10, 12$) over lags 1--60 days, showing strong qualitative stability. Panel B: Real synergy curve for Core H\(\beta\) compared with the median and 95\% percentile envelopes for circular-shift (blue) and block-bootstrap (red) surrogate runs, indicating local significance relative to standard nulls.}
\label{fig:robustness_and_nulls}
\end{figure*}"""

old_fig4 = r"""\begin{figure}[ht]
\centering
\includegraphics[width=\textwidth]{figure4_iccf_vs_surd.png}
\caption{Comparison of the Inter-Correlation Function (ICCF, blue) and the SURD synergy curve (purple) for the Blue Wing (top), Core (middle), and Red Wing (bottom) targets. Vertical lines mark the respective peak lags.}
\label{fig:iccf_vs_surd}
\end{figure}"""

new_fig4 = r"""\begin{figure*}[t]
\centering
\includegraphics[width=0.95\textwidth]{figure4_iccf_vs_surd.png}
\caption{Comparison of the Inter-Correlation Function (ICCF, blue) and the SURD synergy curve (purple) for the Blue Wing (top), Core (middle), and Red Wing (bottom) targets. Vertical lines mark the respective peak lags.}
\label{fig:iccf_vs_surd}
\end{figure*}"""

old_fig5 = r"""\begin{figure}[ht]
\centering
\includegraphics[width=0.8\textwidth]{figure6_javelin_posteriors.png}
\caption{JAVELIN MCMC lag posterior distributions for the Blue Wing (top), Core (middle), and Red Wing (bottom) of H\(\beta\) in NGC~5548. Black solid and dashed lines represent the median and $1\sigma$ (16th--84th percentiles) confidence intervals, respectively.}
\label{fig:javelin_posteriors}
\end{figure}"""

new_fig5 = r"""\begin{figure*}[t]
\centering
\includegraphics[width=0.8\textwidth]{figure6_javelin_posteriors.png}
\caption{JAVELIN MCMC lag posterior distributions for the Blue Wing (top), Core (middle), and Red Wing (bottom) of H\(\beta\) in NGC~5548. Black solid and dashed lines represent the median and $1\sigma$ (16th--84th percentiles) confidence intervals, respectively.}
\label{fig:javelin_posteriors}
\end{figure*}"""

old_fig6 = r"""\begin{figure}[ht]
\centering
\includegraphics[width=\textwidth]{figure5_synthetic_benchmarks.png}
\caption{SURD information decomposition scans for the four synthetic benchmark cases run under realistic NGC~5548 observational conditions (real sampling, seasonal gaps, and flux errors). Curves represent the Monte Carlo median Synergy (purple, with $1\sigma$ shaded envelope), Redundancy (gray dashed), and Unique information for drivers S1 (blue dotted) and S2 (red dash-dot).}
\label{fig:synthetic_benchmarks}
\end{figure}"""

new_fig6 = r"""\begin{figure*}[t]
\centering
\includegraphics[width=0.95\textwidth]{figure5_synthetic_benchmarks.png}
\caption{SURD information decomposition scans for the four synthetic benchmark cases run under realistic NGC~5548 observational conditions (real sampling, seasonal gaps, and flux errors). Curves represent the Monte Carlo median Synergy (purple, with $1\sigma$ shaded envelope), Redundancy (gray dashed), and Unique information for drivers S1 (blue dotted) and S2 (red dash-dot).}
\label{fig:synthetic_benchmarks}
\end{figure*}"""

old_fig7 = r"""\begin{figure}[ht]
\centering
\includegraphics[width=\textwidth]{figure7_history_conditioning.png}
\caption{Comparison of the unconditioned (blue) and target-history conditioned (red) SURD synergy (left) and information leak (right) curves for the Core H\(\beta\) target. Conditioning on the target's own past shifts the synergy peak from a short autocorrelation-dominated delay of 9 days to a robust, long-term cross-component delay of 73 days.}
\label{fig:history_conditioning}
\end{figure}"""

new_fig7 = r"""\begin{figure*}[t]
\centering
\includegraphics[width=0.95\textwidth]{figure7_history_conditioning.png}
\caption{Comparison of the unconditioned (blue) and target-history conditioned (red) SURD synergy (left) and information leak (right) curves for the Core H\(\beta\) target. Conditioning on the target's own past shifts the synergy peak from a short autocorrelation-dominated delay of 9 days to a robust, long-term cross-component delay of 73 days.}
\label{fig:history_conditioning}
\end{figure*}"""

# Apply replacements
text = text.replace(old_fig1, new_fig1)
text = text.replace(old_fig2, new_fig2)
text = text.replace(old_fig3, new_fig3)
text = text.replace(old_fig4, new_fig4)
text = text.replace(old_fig5, new_fig5)
text = text.replace(old_fig6, new_fig6)
text = text.replace(old_fig7, new_fig7)

# Overwrite V3.tex
with open('/Users/ayan/Programs/SURD/overleaf_draft/V3.tex', 'w') as f:
    f.write(text)

print("Finished fixing figures to use figure* and [t] for MNRAS")
