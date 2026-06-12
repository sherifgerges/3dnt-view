# Bundled example datasets

The app shows a **"Try ATP2B2 example"** button only when both of these files exist:

- `ATP2B2_cases.txt`
- `ATP2B2_controls.txt`

Each is a plain list of missense variants, one per line, in either notation
(`E457K` or `p.Glu457Lys`). One line = one allele (or use a tab/comma-delimited
file with an integer allele-count column).

Drop the real ATP2B2 case/control variant lists here with exactly those two
filenames and the example button will appear and run against the cached
AlphaFold model for Q01814.
