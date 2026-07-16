import re

# Load the current V2.tex
with open('/Users/ayan/Programs/SURD/overleaf_draft/V2.tex', 'r') as f:
    content = f.read()

# Replace \textcite with \citet (natbib equivalent)
content = content.replace(r'\textcite', r'\citet')

# Define the new MNRAS preamble and author block
mnras_preamble = r"""\documentclass[fleqn,usenatbib]{mnras}

\usepackage{newtxtext,newtxmath}
\usepackage[T1]{fontenc}
\usepackage{ae,aecompl}

\usepackage{graphicx}
\usepackage{amsmath}

\title[SURD in NGC 5548]{Multivariate and Velocity-Resolved Reverberation in NGC 5548: Astrophysical Motivation for an Information-Based Analysis of Continuum--BLR Coupling}

\author[A. Mitra \& V. Zarikas]{
Ayan Mitra,$^{1}$\thanks{E-mail: ayan@illinois.edu}
Vasilios Zarikas$^{2}$
\\
$^{1}$Department of Astronomy, University of Illinois at Urbana-Champaign, Urbana, IL 61801, USA\\
$^{2}$Department of Physics, University of Thessaly, Lamia 35100, Greece
}

\date{Accepted XXX. Received YYY; in original form ZZZ}
\pubyear{2026}

\begin{document}
\label{firstpage}
\pagerange{\pageref{firstpage}--\pageref{lastpage}}
\maketitle
"""

# Find where the abstract ends
abstract_match = re.search(r'\\begin\{abstract\}(.*?)\\end\{abstract\}', content, re.DOTALL)
if not abstract_match:
    raise ValueError("Could not find abstract in V2.tex")

abstract_text = abstract_match.group(1).strip()

rest_of_doc = content[abstract_match.end():]

# In the rest of the document, we should replace \printbibliography with the BibTeX commands
rest_of_doc = rest_of_doc.replace(r'\printbibliography', r"""\bibliographystyle{mnras}
\bibliography{refs}""")

# Let's check if there is a duplicate \maketitle or \begin{document} in rest_of_doc
rest_of_doc = re.sub(r'\\maketitle', '', rest_of_doc)

new_tex_content = mnras_preamble + "\n\\begin{abstract}\n" + abstract_text + "\n\\end{abstract}\n" + r"""\begin{keywords}
galaxies: active -- galaxies: nuclei -- galaxies: individual: NGC 5548 -- methods: statistical -- methods: observational
\end{keywords}
""" + rest_of_doc

# Add \label{lastpage} before \end{document}
new_tex_content = new_tex_content.replace(r'\end{document}', r"""\label{lastpage}
\end{document}""")

# Overwrite V2.tex directly
with open('/Users/ayan/Programs/SURD/overleaf_draft/V2.tex', 'w') as f:
    f.write(new_tex_content)

print("Conversion script finished overwriting V2.tex")
