with open('/Users/ayan/Programs/SURD/produce_paper_plots.py', 'r') as f:
    plot_text = f.read()

# Change the title in produce_paper_plots.py
plot_text = plot_text.replace(
    "ax.set_title('Figure 8: Spurious Synergy Peaks from Seasonal Windowing (Negative Control)')",
    "ax.set_title('Figure 7: Spurious Synergy Peaks from Seasonal Windowing (Negative Control)')"
)

with open('/Users/ayan/Programs/SURD/produce_paper_plots.py', 'w') as f:
    f.write(plot_text)

print("produce_paper_plots.py updated successfully!")

with open('/Users/ayan/Programs/SURD/overleaf_draft/V3.tex', 'r') as f:
    v3_text = f.read()

# Change the contradiction in V3.tex
old_phrase = "while the synergy peaks at 56d, 109d, and 159d locally exceed the surrogate envelopes in uncorrected scans, they do not survive global multiple-testing controls."
new_phrase = "while the synergy peaks at 56d and 109d locally exceed the surrogate envelopes in uncorrected scans (and the 159d Core peak does not clear even local surrogate thresholds), none of these peaks survive global multiple-testing controls."

v3_text = v3_text.replace(old_phrase, new_phrase)

with open('/Users/ayan/Programs/SURD/overleaf_draft/V3.tex', 'w') as f:
    f.write(v3_text)

print("V3.tex updated successfully!")
