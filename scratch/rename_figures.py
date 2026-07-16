import os

# 1. Rename files in overleaf_draft
draft_dir = '/Users/ayan/Programs/SURD/overleaf_draft/'
h_cond_old = os.path.join(draft_dir, 'figure7_history_conditioning.png')
h_cond_new = os.path.join(draft_dir, 'figure8_history_conditioning.png')
s_alias_old = os.path.join(draft_dir, 'figure8_seasonal_aliasing.png')
s_alias_new = os.path.join(draft_dir, 'figure7_seasonal_aliasing.png')

# Swap them carefully
temp_h_cond = os.path.join(draft_dir, 'temp_h_cond.png')
if os.path.exists(h_cond_old):
    os.rename(h_cond_old, temp_h_cond)
if os.path.exists(s_alias_old):
    os.rename(s_alias_old, s_alias_new)
if os.path.exists(temp_h_cond):
    os.rename(temp_h_cond, h_cond_new)

print("Renamed files in overleaf_draft.")

# 2. Update V3.tex
with open('/Users/ayan/Programs/SURD/overleaf_draft/V3.tex', 'r') as f:
    v3_text = f.read()

# Swap the image files
v3_text = v3_text.replace('figure7_history_conditioning.png', 'figure8_history_conditioning.png')
v3_text = v3_text.replace('figure8_seasonal_aliasing.png', 'figure7_seasonal_aliasing.png')

with open('/Users/ayan/Programs/SURD/overleaf_draft/V3.tex', 'w') as f:
    f.write(v3_text)

print("Updated V3.tex with new image filenames.")

# 3. Update produce_paper_plots.py
with open('/Users/ayan/Programs/SURD/produce_paper_plots.py', 'r') as f:
    plot_text = f.read()

# Swap in plot script
plot_text = plot_text.replace('figure7_history_conditioning.png', 'figure8_history_conditioning.png')
plot_text = plot_text.replace('figure8_seasonal_aliasing.png', 'figure7_seasonal_aliasing.png')

# Also swap the labels in comments/prints if any
plot_text = plot_text.replace("Generating Figure 7: Target-History Conditioning", "Generating Figure 8: Target-History Conditioning")
plot_text = plot_text.replace("Figure 7: figure7_history_conditioning.png", "Figure 8: figure8_history_conditioning.png")
plot_text = plot_text.replace("Generating Figure 8: Seasonal Aliasing", "Generating Figure 7: Seasonal Aliasing")
plot_text = plot_text.replace("Figure 8: figure8_seasonal_aliasing.png", "Figure 7: figure7_seasonal_aliasing.png")

with open('/Users/ayan/Programs/SURD/produce_paper_plots.py', 'w') as f:
    f.write(plot_text)

print("Updated produce_paper_plots.py with new image filenames.")
