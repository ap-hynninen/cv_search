#!/bin/bash

grep --no-filename -f residues.txt step1_first10k_pdbs/1.run/???.* > coord1.txt
grep --no-filename -f residues.txt step1_first10k_pdbs/1.run/????.* >> coord1.txt
grep --no-filename -f residues.txt step1_first10k_pdbs/2.run/???.* > coord2.txt
grep --no-filename -f residues.txt step1_first10k_pdbs/2.run/????.* >> coord2.txt
grep --no-filename -f residues.txt step1_first10k_pdbs/3.run/???.* > coord3.txt
grep --no-filename -f residues.txt step1_first10k_pdbs/3.run/????.* >> coord3.txt
grep --no-filename -f residues.txt step1_first10k_pdbs/4.run/???.* > coord4.txt
grep --no-filename -f residues.txt step1_first10k_pdbs/4.run/????.* >> coord4.txt
grep --no-filename -f residues.txt step1_first10k_pdbs/5.run/???.* > coord5.txt
grep --no-filename -f residues.txt step1_first10k_pdbs/5.run/????.* >> coord5.txt

grep --no-filename -f residues.txt step1_second10k_pdbs/1/???.* > coord6.txt
grep --no-filename -f residues.txt step1_second10k_pdbs/1/????.* >> coord6.txt
grep --no-filename -f residues.txt step1_second10k_pdbs/2/???.* > coord7.txt
grep --no-filename -f residues.txt step1_second10k_pdbs/2/????.* >> coord7.txt
grep --no-filename -f residues.txt step1_second10k_pdbs/3/???.* > coord8.txt
grep --no-filename -f residues.txt step1_second10k_pdbs/3/????.* >> coord8.txt
grep --no-filename -f residues.txt step1_second10k_pdbs/4/???.* > coord9.txt
grep --no-filename -f residues.txt step1_second10k_pdbs/4/????.* >> coord9.txt
grep --no-filename -f residues.txt step1_second10k_pdbs/5/???.* > coord10.txt
grep --no-filename -f residues.txt step1_second10k_pdbs/5/????.* >> coord10.txt

cat coord1.txt > coords.txt
cat coord2.txt >> coords.txt
cat coord3.txt >> coords.txt
cat coord4.txt >> coords.txt
cat coord5.txt >> coords.txt

cat coord6.txt >> coords.txt
cat coord7.txt >> coords.txt
cat coord8.txt >> coords.txt
cat coord9.txt >> coords.txt
cat coord10.txt >> coords.txt
