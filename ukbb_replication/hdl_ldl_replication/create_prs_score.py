# Given PLINK genotypes for a trait, construct a PolyGenic Risk (PGR) Score

import argparse


def parseargs():    # handle user arguments
	parser = argparse.ArgumentParser(description="Given PLINK genotypes for " +
		"a trait, construct a PolyGenic Risk (PGR) Score.")
	parser.add_argument('--ped', required = True, help = 'Plink PED file.')
	parser.add_argument('--map', required = True, help = 'Plink MAP file.')
	parser.add_argument('--effect_alleles', required = True,
		help = 'tsv file with columns: rsid, effect-allele, beta.')
	parser.add_argument('--output', required = True, help = 'Where to write PGR scores.')
	args = parser.parse_args()
	return args


def read_map(mapfile):
	ordered_snps = []
	with(open(mapfile, 'r')) as infile:
		for line in infile:
			ordered_snps.append(line.split('\t')[1])  # rsids
	return ordered_snps


def read_effect_alleles(effect_alleles_file):
	snp_to_effect_allele = {}
	with(open(effect_alleles_file, 'r')) as infile:
		for line in infile:
			splits = line.strip().split('\t')
			rsid, effect_allele, beta = splits
			beta = float(beta)
			snp_to_effect_allele[rsid] = [effect_allele, beta]
	return snp_to_effect_allele


def calc_pgr_score(genos, ordered_snps, snp_to_effect_allele):
	pgr = 0.0
	for i in range(len(ordered_snps)):
		snp_loc = ordered_snps[i]
		copy1, copy2 = genos[i * 2], genos[i * 2 + 1]
		effect_allele, beta = snp_to_effect_allele[snp_loc]
		if copy1 == effect_allele:
			pgr += beta
		if copy2 == effect_allele:
			pgr += beta
	return str(pgr)


def calculate_scores(pedfile, outputfile, ordered_snps, snp_to_effect_allele):
	with(open(pedfile, 'r')) as infile:
		with(open(outputfile, 'w')) as outfile:
			if line.startswith('#'):
				continue
			for line in infile:
				splits = line.strip().split()
				iid, genos = splits[1], splits[6:]
				pgr = calc_pgr_score(genos, ordered_snps, snp_to_effect_allele)
				outfile.write(iid + '\t' + pgr + '\n')


def main():
	args = parseargs()
	ordered_snps = read_map(args.map)
	snp_to_effect_allele = read_effect_alleles(args.effect_alleles)
	calculate_scores(args.ped, args.output, ordered_snps, snp_to_effect_allele)


if __name__ == '__main__':
	main()
#

