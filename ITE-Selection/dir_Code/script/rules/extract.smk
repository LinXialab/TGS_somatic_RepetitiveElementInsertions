rule mask_trf:
    input:
        get_full_path
    output:
        config["path"]["output_trf"] + "/{ins_id}.out"
    params:
        tool=config["tool"]["trf"],
        trf=config["parameter"]["trf"],
    shell:
        "{params.tool} {input} {params.trf} > {output}"

rule mask_dust:
    input:
        get_full_path
    output:
        config["path"]["output_sdust"] + "/{ins_id}.out"
    params:
        config["tool"]["sdust"]
    shell:
        "{params} {input} > {output}"

rule summary_repeatmasker:
    output:
        config["file"]["repeatmasker_summary"]
    run:
        cor_dict = dict()
        type_cor_dict = dict()
        result_dict = dict()
        for ins_id in repeat_ls:
            filepath = get_repeat_path(ins_id)
            out_df = parse_repeat(filepath)
            if isinstance(out_df, int):
                result_dict.setdefault(ins_id, 0)
                cor_dict.setdefault(ins_id, 0)
                type_cor_dict.setdefault(ins_id, 0)
            else:
                family_ls = out_df["TEfamily"].tolist()
                len_ls = out_df["ilength"].tolist()
                len_s = "|".join(["{}_{}".format(len_ls[i], family_ls[i]) for i in range(len(family_ls))])
                result_dict.setdefault(ins_id, len_s)
                cor_dict.setdefault(ins_id, "|".join(out_df["family_cor"].tolist()))
                type_cor_dict.setdefault(ins_id, "|".join(out_df["type_cor"].tolist()))
        pd.DataFrame({"len": result_dict, "coordinate": cor_dict, "type_cor": type_cor_dict}).to_csv(output[0], sep="\t") 

rule check_trf_fasta:
    input:
        trf=expand(config["path"]["output_trf"] + "/{ins_id}.out", ins_id=iter(sv_ls)),
        sdust=expand(config["path"]["output_sdust"] + "/{ins_id}.out", ins_id=iter(sv_ls))
    output:
        config["file"]["trf_size_check"]
    run:
        os.mkdir(config["path"]["pattern_f"])
        valid_sv = []
        for ins_id in sv_ls:
            trf_ls = trf_seq("{}/{}.out".format(config["path"]["output_trf"], ins_id))
            if not trf_ls:
                continue
            valid_sv.append(ins_id)
            output_p = "{}/{}.fasta".format(config["path"]["pattern_f"], ins_id)
            with open(output_p, "w+") as handle:
                for idx, seq_region in enumerate(trf_ls):
                    SeqIO.write(SeqRecord(Seq(seq_region[-1]), id="{}_{}".format(ins_id, idx), description=""), handle, "fasta")
        pd.Series(valid_sv).to_csv(output[0], index=False, header=False)

