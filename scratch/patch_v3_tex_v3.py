with open('/Users/ayan/Programs/SURD/overleaf_draft/V3.tex', 'r') as f:
    text = f.read()

# Replacement 1: First sentence of Section 4.1
old_1 = "To test the physical stability of the synergy and leak profiles, we extended the lag search range from the initial limit of 60 days to 120 days."
new_1 = "To test the physical stability of the synergy and leak profiles, we extended the lag search range from the initial limit of 60 days to 200 days."
text = text.replace(old_1, new_1)

# Replacement 2: Section 4.2 surrogate test numbers and FDR bolding
old_2 = r"""At the peak synergy lags identified in the extended scan:
\begin{itemize}
\item \textbf{Core H\(\beta\) (at 159 days):} Real unnormalized synergy is $0.0526$, compared to the 95th percentile surrogate thresholds from the circular shift ($\sim 0.136$) and block bootstrap ($\sim 0.086$) null distributions.
\item \textbf{Red Wing H\(\beta\) (at 56 days):} Real unnormalized synergy is $0.5677$.
\item \textbf{Blue Wing H\(\beta\) (at 109 days):} Real unnormalized synergy is $0.2198$.
\end{itemize}

While these real synergy values locally exceed the 95th percentile null envelopes at specific lags, a key statistical caveat must be emphasized: **after applying the Benjamini–Hochberg False Discovery Rate (FDR) correction for multiple testing across the 200 scanned lags, none of the individual lags remain globally statistically significant.**"""

new_2 = r"""At the peak synergy lags identified in the extended scan, the real unnormalized synergy values are $0.0526$~bits for the Core (at 159 days), $0.5677$~bits for the Red Wing (at 56 days), and $0.2198$~bits for the Blue Wing (at 109 days). For the standard lag range (1--60 days), the real unnormalized synergy curves for all three targets locally exceed the 95th percentile surrogate null envelopes (for example, the Core unnormalized synergy at the standard 45-day peak is $0.1883$~bits, which is well above the circular shift threshold of $\sim 0.136$~bits and the block bootstrap threshold of $\sim 0.086$~bits). However, in the extended scan (up to 200 days), while the Red and Blue wing unnormalized peaks exceed the local surrogate null envelopes, the Core peak at 159 days does not.

Crucially, after applying the Benjamini–Hochberg False Discovery Rate (FDR) correction for multiple testing across the 200 scanned lags, \textbf{none of the individual lags remain globally statistically significant for any of the three targets.}"""

text = text.replace(old_2, new_2)

# Replacement 3: Section 4.4 closing sentence
old_3 = "In contrast, the much longer-term SURD synergy peaks (76--119 days) suggest the presence of a broader multivariate information structure that is invisible to these standard linear frameworks."
new_3 = "In contrast, the long-lag unconditioned SURD synergy peaks (56--159 days) are consistent with windowing and seasonal aliasing artifacts (as demonstrated in Section 4.7) rather than physical BLR structures."
text = text.replace(old_3, new_3)

# Replacement 4: Section 4.5 FDR and global maximum-statistic section
old_4 = r"""Additionally, we applied a False Discovery Rate (FDR) multiple-testing correction to the SURD lag scans. Because we scan $M=120$ individual lags across multiple components, the probability of obtaining false-positive significance peaks increases. We estimated the p-value at each lag by comparing the real synergy to the normal-approximated cont-shuffle null distribution:
\begin{equation}
p(\tau) = 1 - \Phi\left(\frac{S_{12}(\tau) - \mu_{\mathrm{null}}}{\sigma_{\mathrm{null}}}\right),
\end{equation}
where $\Phi$ is the standard normal CDF, and $\mu_{\mathrm{null}}$ and $\sigma_{\mathrm{null}}$ are the median and standard deviation of the surrogate distribution. Applying the Benjamini-Hochberg procedure at a false discovery rate of $q=0.05$ reveals that \textbf{none of the individual lags are statistically significant under a global multiple-comparison control}.

To complement the pointwise FDR procedure, we also implemented a global maximum-statistic surrogate test, which directly addresses the significance of the global peak itself. For each of $N_{\mathrm{surrogate}} = 500$ circular shift surrogate realizations, we compute the SURD scan over the full range of 120 lags and record the maximum synergy value ($S_{\mathrm{max, null}}$). This yields an empirical distribution of the peak synergy values expected purely by chance under the global null hypothesis. The empirical global p-value is calculated as:
\begin{equation}
p_{\mathrm{global}} = \frac{1 + \#(S_{\mathrm{max, null}} \ge S_{\mathrm{max, real}})}{1 + N_{\mathrm{surrogate}}}.
\end{equation}
Comparing our real maximum synergy values ($0.2284$ for Core, $0.3386$ for Red, and $0.2843$ for Blue) to this global max-statistic null distribution yields global p-values of $p_{\mathrm{global}} = 0.9960$ (Core), $p_{\mathrm{global}} = 0.6547$ (Red), and $p_{\mathrm{global}} = 0.7725$ (Blue). This global test indicates that the observed synergy peaks are consistent with random fluctuations of scanned red-noise curves under a global null procedure, supporting the cautionary conclusion that the synergy profiles represent candidate coupling envelopes requiring validation in independent data rather than precise physical delays.

This is a major methodological contribution for the paper: while the synergy peaks at 76d, 81d, and 119d stand out as local features that locally exceed the surrogate envelopes ($p < 0.05$ individually), they must be interpreted as broad lag-dependent information profiles (or candidate coupling envelopes requiring validation in independent data) rather than precise, delta-function-like delay times. Multiple-testing corrections are essential for information-theoretic scans of variable astronomical light curves to avoid over-interpreting isolated lag values."""

new_4 = r"""Additionally, we applied a False Discovery Rate (FDR) multiple-testing correction to the SURD lag scans. Because we scan $M=200$ individual lags across multiple components, the probability of obtaining false-positive significance peaks increases. We estimated the p-value at each lag by comparing the real synergy to the normal-approximated cont-shuffle null distribution:
\begin{equation}
p(\tau) = 1 - \Phi\left(\frac{S_{12}(\tau) - \mu_{\mathrm{null}}}{\sigma_{\mathrm{null}}}\right),
\end{equation}
where $\Phi$ is the standard normal CDF, and $\mu_{\mathrm{null}}$ and $\sigma_{\mathrm{null}}$ are the median and standard deviation of the surrogate distribution. Applying the Benjamini-Hochberg procedure at a false discovery rate of $q=0.05$ reveals that \textbf{none of the individual lags are statistically significant under a global multiple-comparison control}.

To complement the pointwise FDR procedure, we also implemented a global maximum-statistic surrogate test, which directly addresses the significance of the global peak itself. For each of $N_{\mathrm{surrogate}} = 100$ circular shift surrogate realizations, we compute the SURD scan over the full range of 200 lags and record the maximum normalized synergy value ($\widehat{S}_{\mathrm{max, null}}$). This yields an empirical distribution of the peak normalized synergy values expected purely by chance under the global null hypothesis. The empirical global p-value is calculated as:
\begin{equation}
p_{\mathrm{global}} = \frac{1 + \#(\widehat{S}_{\mathrm{max, null}} \ge \widehat{S}_{\mathrm{max, real}})}{1 + N_{\mathrm{surrogate}}}.
\end{equation}
Comparing our real maximum normalized synergy values ($0.5127$ for Core, $0.6390$ for Red, and $0.5312$ for Blue) to this global max-statistic null distribution yields global p-values of $p_{\mathrm{global}} = 1.00$ (Core), $p_{\mathrm{global}} = 0.63$ (Red), and $p_{\mathrm{global}} = 0.98$ (Blue). This global test indicates that the observed synergy peaks are consistent with random fluctuations of scanned red-noise curves under a global null procedure, supporting the conclusion that the unconditioned synergy peaks are non-significant windowing and aliasing artifacts rather than real delays.

This is a major methodological contribution for the paper: while the synergy peaks at 56d, 109d, and 159d locally exceed the surrogate envelopes in uncorrected scans, they do not survive global multiple-testing controls. Multiple-testing corrections and negative controls are essential for information-theoretic scans of variable astronomical light curves to avoid over-interpreting isolated lag values."""

text = text.replace(old_4, new_4)

# Replacement 5: Section 5.1 opening (stale numbers and co-equal hypotheses)
old_5 = r"""The results in Section~\ref{sec:results} reveal a notable discrepancy between classical linear delays (13--20 days, as measured by both ICCF and JAVELIN) and the SURD synergy peaks (76--119 days). In this section, we discuss the physical implications of these findings for the Broad-Line Region (BLR) of NGC~5548.

\subsection{The Multi-Component Broad-Line Region}

The H\(\beta\) reverberation lag of NGC~5548 is historically reported in the range of 10 to 20 days (e.g., $15.6 \pm 0.5$ days during the high-cadence AGN STORM campaign \cite{Pei2017}). Our ICCF peaks of 13.0 days (Red Wing), 19.0 days (Blue Wing), and 20.0 days (Core), as well as JAVELIN peak delays of $12.8 \pm 1.8$ days (Red), $19.4_{-2.9}^{+2.5}$ days (Blue), and $20.1_{-3.2}^{+2.9}$ days (Core), are in excellent agreement with these standard measurements. They reflect the prompt, linear light-travel response of the broad-line gas to continuum fluctuations.

The fact that the SURD synergy peaks occur at much longer timescales (76 days for the Core, 81 days for the Blue Wing, and 119 days for the Red Wing) suggests that synergy is sensitive to different physical mechanisms or systematic effects than simple linear reprocessing:
\begin{enumerate}
\item \textbf{Extended Geometry:} The long-lag synergy could trace the outer, slow-moving regions of the BLR. In an extended, disk-like or virialized BLR, the outer regions respond with a longer light-travel delay and contribute to the joint information structure of the line profile.
\item \textbf{Non-linear Reprocessing:} While the ICCF measures linear correlation, SURD is sensitive to non-linear statistical dependencies. The H\(\beta\) emissivity is known to depend non-linearly on the ionizing flux and local gas density, which can manifest as delayed synergistic information that requires both the driving continuum and other line components to resolve.
\item \textbf{Red Noise and Seasonal Aliasing:} Alternatively, as demonstrated by our synthetic negative control test, the combination of seasonal observing gaps, linear interpolation (which smooths high-frequency variance), and normalized synergy inflation as joint mutual information approaches zero artificially manufactures spurious peaks at long lags. These peaks do not pass global significance tests (FDR) and are scan-range-dependent, meaning they represent windowing and boundary artifacts rather than physical coupling.
\end{enumerate}"""

new_5 = r"""The results in Section~\ref{sec:results} reveal a notable discrepancy between classical linear delays (13--20 days, as measured by both ICCF and JAVELIN) and the SURD synergy peaks (56--159 days). In this section, we discuss the physical implications of these findings for the Broad-Line Region (BLR) of NGC~5548.

\subsection{The Multi-Component Broad-Line Region}

The H\(\beta\) reverberation lag of NGC~5548 is historically reported in the range of 10 to 20 days (e.g., $15.6 \pm 0.5$ days during the high-cadence AGN STORM campaign \cite{Pei2017}). Our ICCF peaks of 13.0 days (Red Wing), 19.0 days (Blue Wing), and 20.0 days (Core), as well as JAVELIN peak delays of $12.8 \pm 1.8$ days (Red), $19.4_{-2.9}^{+2.5}$ days (Blue), and $20.1_{-3.2}^{+2.9}$ days (Core), are in excellent agreement with these standard measurements. They reflect the prompt, linear light-travel response of the broad-line gas to continuum fluctuations.

The fact that the SURD synergy peaks occur at much longer timescales (159 days for the Core, 109 days for the Blue Wing, and 56 days for the Red Wing) suggests that these unconditioned peaks must be interpreted against systematic windowing and aliasing effects rather than as direct physical delays:
\begin{enumerate}
\item \textbf{Extended Geometry:} While extended geometry (tracing the outer regions of the BLR with longer light-travel delays) or non-linear photoionization reprocessing could theoretically contribute to long-lag information coupling, our synthetic negative control (Section 4.7) demonstrates that seasonal windowing, interpolation-induced red-noise coherence, and normalized synergy inflation are fully sufficient to manufacture the observed peaks without invoking any physical long-lag structure.
\item \textbf{Non-linear Reprocessing:} Similarly, non-linear reprocessing is rendered redundant as a necessary explanation, since the spurious peaks in our control simulations achieve amplitudes ($\widehat{S}_{12} \sim 0.59$) that fully match or exceed those in the real data.
\item \textbf{Red Noise and Seasonal Aliasing:} Consequently, the combination of seasonal observing gaps, linear interpolation, and normalized synergy inflation is the primary, globally consistent explanation for the unconditioned peaks. This explains why they are highly scan-range-dependent and fail to survive global multiple-testing controls.
\end{enumerate}"""

text = text.replace(old_5, new_5)

# Replacement 6: Section 5.2 closing sentence
old_6 = "As discussed in Section~\\ref{sec:limitations}, this mixing of physical regimes can produce artificial state-identification synergy and explains why the broad 76--119 day synergy envelopes do not survive global multiple-testing controls."
new_6 = "As discussed in Section~\\ref{sec:limitations}, this mixing of physical regimes may compound the windowing and aliasing artifacts, although the negative control in Section 4.7 alone is sufficient to explain the observed peak amplitudes and why they do not survive global multiple-testing controls."
text = text.replace(old_6, new_6)

# Replacement 7: Escape % in conclusions
old_7 = "check (89.2% empty cells) reveal"
new_7 = r"check (89.2\% empty cells) reveal"
text = text.replace(old_7, new_7)

with open('/Users/ayan/Programs/SURD/overleaf_draft/V3.tex', 'w') as f:
    f.write(text)

print("V3.tex patched successfully!")
