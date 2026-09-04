import gwaslab as gl 

ss = gl.Sumstats("mama_gp2Only_non-eur.gwaslab.tsv.gz", build="hg38", fmt="gwaslab")

ss.fix_chr(remove=True)
ss.fix_pos(remove=True)

res = ss.get_novel(
            known="./known_lociPD.tsv",
            only_novel=False,            # return both
            sig_level=5e-8,
            windowsizekb=500,
            windowsizekb_for_novel=500,
anno=True
        )

res.query("NOVEL == True").to_csv("out_file_novel", sep="\t", index=False)
res.query("NOVEL == False").to_csv("out_file_known", sep="\t", index=False)
