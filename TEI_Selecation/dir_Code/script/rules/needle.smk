rule extract_needle:
    input:
        get_full_path
    output:
        dw=config["path"]["ref_f"] + "/{ins_id}_dw.fasta",
        up=config["path"]["ref_f"] + "/{ins_id}_up.fasta"
    params:
        "{ins_id}"
    run:
        records = SeqIO.parse(input[0], "fasta")
        record = list(records)[0]
        total_len = len(record.seq)
        ins_site = insert_df.loc[params[0], "ins_loc_hg38"]
        ref_info = ins_site.split("_")
        ref_chr, ref_start = ref_info[0], int(ref_info[1])
        ref_end = ref_start+1
        #ref_length = {}
        #for seq_record in ref_records:
        #    ref_length[str(seq_record.id)] = len(seq_record)
	#ref_length
        if ref_start+1-(total_len+50) < 0:
           new_up_start = 0
        else:
           new_up_start = ref_start+1-(total_len+50)
        if ref_end+total_len+50 > len(ref_records_dict[ref_chr].seq):
           new_dw_end = len(ref_records_dict[ref_chr].seq)
        else:
           new_dw_end = ref_end+total_len+50
        ref_up_seq = ref_records_dict[ref_chr].seq[new_up_start:ref_start+1]
        ref_dw_seq = ref_records_dict[ref_chr].seq[ref_end:new_dw_end]
        with open(output.up, "w+") as handle:
            SeqIO.write(SeqRecord(ref_up_seq, id="{}_up".format(params[0]), description=""), handle, "fasta")
        with open(output.dw, "w+") as handle:
            SeqIO.write(SeqRecord(ref_dw_seq, id="{}_dw".format(params[0]), description=""), handle, "fasta")

rule needle_up:
    input:
        q=get_full_path,
        s=config["path"]["ref_f"] + "/{ins_id}_up.fasta"
    output:
        config["path"]["output_needle"] + "/{ins_id}_up.out"
    params:
        tool=config["tool"]["needle"],
        needle=config["parameter"]["needle"]
    shell:
        "{params.tool} -asequence {input.q} -bsequence {input.s} {params.needle} -out {output}"

rule needle_dw:
    input:
        q=get_full_path,
        s=config["path"]["ref_f"] + "/{ins_id}_dw.fasta"
    output:
        config["path"]["output_needle"] + "/{ins_id}_dw.out"
    params:
        tool=config["tool"]["needle"],
        needle=config["parameter"]["needle"]
    shell:
        "{params.tool} -asequence {input.q} -bsequence {input.s} {params.needle} -out {output}"

rule summary_needle:
    input:
        expand(config["path"]["output_needle"] + "/{ins_id}_{direc}.out", ins_id=iter(sv_ls), direc=iter(["up", "dw"]))
    output:
        config["file"]["needle_summary"]
    run:
        result_dict = dict()
        needle_f = config["path"]["output_needle"]
        for filename in os.listdir(needle_f):
            filepath = os.path.join(needle_f, filename)
            region = filename.split(".")[-2].split("_")[-1]
            ins_id = "_".join(filename.split(".")[0].split("_")[:-1])
            
            align = AlignIO.read(filepath, "emboss")
            q_ungap_start, q_ungap_end = get_ungap_cor(align[0])
            s_ungap_start, s_ungap_end = get_ungap_cor(align[1])
            align_len = len(set(range(q_ungap_start, q_ungap_end)).intersection(range(s_ungap_start, s_ungap_end)))
            q_len = len(str(align[0].seq.ungap("-")))
            align_q_cov = align_len/q_len
            iden_perc = align_len and align.annotations["identity"]/align_len or 0
            result_dict.setdefault(ins_id, {}).update({"{}_align_q_cov".format(region): align_q_cov,
                                                       "{}_iden_perc".format(region): iden_perc,
                                                       "{}_start_dis".format(region): q_ungap_start-s_ungap_start,
                                                       "{}_end_dis".format(region): q_ungap_end-s_ungap_end})

        pd.DataFrame(result_dict).T.to_csv(output[0], sep="\t")
