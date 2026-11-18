#!/usr/bin/env python
"""Compare the numbers that the manuscript text expects with the numbers that the analysis produced.

    usage:  python checks/check_expected_values.py <manuscript folder>            report differences
            python checks/check_expected_values.py <manuscript folder> --update   rewrite the comment blocks

Every .tex file in <manuscript folder>/bioapp/ starts with a comment block

    % --- expected values (written by checks/check_expected_values.py --update) ---
    %   \\FRFmhDnds = 1.85      FracMinHash dN/dS of the pair, all genes
    % --- end expected values ---

that lists each number macro used in the file and the value it had when the text was written.
The values are compared with the macros in output/numbers/ of this repository. After re-running
the analysis, this script reports which values have changed, so that the sentences that quote
them can be re-read before the new numbers are exported to the manuscript (export_to_manuscript.sh).
It also reports macros that are used in the text but not defined.
"""
import sys
import re
import glob

BEGIN = "% --- expected values (written by checks/check_expected_values.py --update) ---"
END = "% --- end expected values ---"
manuscript = sys.argv[1].rstrip("/")
update = "--update" in sys.argv

# every macro defined in output/numbers/:  name -> (value, meaning)
defined = {}
for path in glob.glob("output/numbers/*.tex"):
    for line in open(path):
        match = re.match(r"\\providecommand\{\\(\w+)\}\{(.*?)\}\s*%\s*(.*)", line)
        if match:
            defined[match.group(1)] = (match.group(2), match.group(3).strip())
prefixes = sorted({name[:2] for name in defined})

problems = 0
for path in sorted(glob.glob(f"{manuscript}/bioapp/*.tex")):
    text = open(path).read()
    if BEGIN not in text:
        continue
    head, rest = text.split(BEGIN, 1)
    block, body = rest.split(END, 1)
    # macros used in the file (also inside comments) that start with one of the number prefixes
    used = sorted(set(re.findall(r"\\((?:%s)[A-Z][A-Za-z]+)" % "|".join(prefixes), body)), key=body.find)
    expected = dict(re.findall(r"%\s+\\(\w+) = (.*?)\s{4,}", block))
    for name in used:
        if name not in defined:
            print(f"{path}: \\{name} is used but not defined in output/numbers/")
            problems += 1
        elif not update and expected.get(name) != defined[name][0]:
            print(f"{path}: \\{name} expected {expected.get(name, '(not listed)')}, now {defined[name][0]}   [{defined[name][1]}]")
            problems += 1
    if update:
        lines = [f"%   \\{name} = {defined[name][0]}      {defined[name][1]}" for name in used if name in defined]
        open(path, "w").write(head + BEGIN + "\n" + "\n".join(lines) + "\n" + END + body)
print("expected-value blocks rewritten" if update else f"{problems} difference(s) between expected and current values")
sys.exit(1 if problems and not update else 0)
