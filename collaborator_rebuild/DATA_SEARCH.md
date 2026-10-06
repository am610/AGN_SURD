# Data search and recovered settings

## Local findings

`history/documents/Tasks.docx` and `history/documents/Tasks_2.docx`
contain an appendix with the following
recommended common lag ranges. These are local analysis recommendations,
not measurements published by Xi.

| Season | Common lag range in days |
| :--- | :--- |
| 2015 | 0 through 20 |
| 2018 | 0 through 20 |
| 2019 | 0 through 20 |
| 2020 | 0 through 25 |
| 2021 | 0 through 20 |
| 2023 | 0 through 30 |

The appendix recommends signed component differential delays from minus 15
through plus 15 days. It discusses fixed candidate bins with blue from minus
3000 to minus 1000, core from minus 1000 to plus 1000, and red from plus 1000
to plus 3000 kilometres per second. The document explicitly asks for a review
of the velocity choice. These are therefore candidate project settings rather
than the authors' exact published bins.

The local AGN Watch spectra span the historical campaign around 1988 through
1993. They cannot supply missing observations in the six modern seasons.

## Public sources checked

Lu Table 2 is available from
https://cdsarc.cds.unistra.fr/ftp/J/ApJS/263/10/table2.dat
with the catalogue description at
https://cdsarc.cds.unistra.fr/ftp/J/ApJS/263/10/ReadMe

Xi Table 1 was retrieved through the machine readable link on the publisher's
article page at https://iopscience.iop.org/article/10.3847/1538-4357/ae1ccb
The actual table URL is
https://content.cld.iop.org/journals/0004-637X/995/2/157/revision1/apjae1ccbt1_mrt.txt

The Lu and Xi article source archives were inspected in the preceding session.
The Lu archive includes its integrated light curve table. The Xi archive
contains article text and figures. Neither inspected archive supplies the
individual broad H beta spectra or blue, core, and red time series.

Targeted searches for the two DOI identifiers, campaign spectra, and public
data repositories did not locate those component data. The guessed CDS
catalogue path for Xi returned HTTP 404. The publisher's dedicated data page
returned a bot verification page, although the article and Table 1 download
succeeded. Accordingly this search does not establish that no public spectra
exist anywhere.

## Paper definitions

Xi section 5.2 describes integration of net broad H beta profiles in a window
from 4700 to 5050 Angstrom and uses 15 and 30 bins of equal flux. Those bins
are different from the three fixed candidate project bins. Published velocity
lag summaries do not determine the underlying flux at every observation.

Lu Table 1 defines campaign date boundaries. Its 2019 season begins on
28 November 2018 and its 2021 season begins on 24 December 2020. The corrected
counts are 62, 40, 81, 52, and 80 for the five Lu seasons. The preceding
calendar year split was incorrect and has been replaced.

The complete Xi campaign has 74 observations and starts in December 2022.
Its campaign label is 2023. Splitting it by Gregorian year would corrupt this
campaign definition. The public preparation script preserves the full table.

## Remaining input gap

F5100 and integrated broad H beta are now prepared for all six seasons, with
389 observations in total. Blue, core, and red H beta remain unavailable in
the sources inspected. Exact predictor specific delays are also absent;
the local appendix explicitly says the papers do not tabulate coarse component
delays. The full requested analysis requires the component curves or calibrated
net broad H beta profiles. These should include uncertainties and wavelength
and flux conventions.

The published author contact for the Xi campaign is xiwenzhe@ynao.ac.cn,
listed by Yunnan Observatories at
https://english.ynao.ac.cn/research/rp/202602/t20260210_1150641.html
No external data request has been sent.
