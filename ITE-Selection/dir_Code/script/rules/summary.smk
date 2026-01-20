rule summary_repeat_check:
    input:
        tf_s=config["file"]["sdust_trf_summary"],
        ne_s=config["file"]["needle_summary"],
    output:
        summary_dict=config["file"]["repeat_check"]["summary_dict"],
        cor_dict=config["file"]["repeat_check"]["cor_dict"],
        ins_ls=config["file"]["repeat_check"]["ins_ls"],
    run:
        result_dict = dict()
        valid_sv_ls = pd.read_csv(config["file"]["trf_size_check"], names=["id"])["id"].tolist()
        for ins_id in valid_sv_ls:
            trf_path = "{}/{}.out".format(config["path"]["output_trf"], ins_id)            
            trf_ls = trf_seq(trf_path)
            ## remove trf annotation with more than one pattern
            if len(trf_ls) > 1:
                continue
            ## remove ins seq trf coverage less than 90%
            fasta_path = get_fasta_path(ins_id)
            records = SeqIO.parse(fasta_path, "fasta")
            record = list(records)[0]
            total_len = len(record.seq)
            trf_len = trf_ls[0][1]-trf_ls[0][0]+1
            if trf_len < total_len*0.9:
                continue
            up_fasta = "{}/{}_up.fasta".format(config["path"]["ref_f"], ins_id)
            up_seq = str(list(SeqIO.parse(up_fasta, "fasta"))[0].seq)
            dw_fasta = "{}/{}_dw.fasta".format(config["path"]["ref_f"], ins_id)
            dw_seq = str(list(SeqIO.parse(dw_fasta, "fasta"))[0].seq)
            pattern_len = len(trf_ls[0][-1])
            extract_len = max([10, int(pattern_len*1.5)])+10
            result_dict.setdefault(ins_id, {}).update({"ins_pattern": trf_ls[0][-1],
                                                       "ins_len": total_len,
                                                       "ins_annot_cor": "_".join(map(str,trf_ls[0][:-1])),
                                                       "up_seq": up_seq[-min([len(up_seq), extract_len]):],
                                                       "dw_seq": dw_seq[:min([len(dw_seq), extract_len])]})
        with open(output.summary_dict, "wb") as handle:
            pickle.dump(result_dict, handle)
        ## get align coordinate
        cor_dict = get_cor_dict_from_p_align(result_dict)
        with open(output.cor_dict, "wb") as handle:
            pickle.dump(cor_dict, handle)
        ## filter alignment by score 
        aligner = Align.PairwiseAligner()
        aligner.query_internal_gap_score = -10
        aligner.query_internal_extend_gap_score = -0.5
        aligner.mismatch_score = -1
        aligner.mode = "local"
        ins_ls = []
        for ins_id, ins_dict in result_dict.items():
            up_similar_score = cal_align_score(ins_dict["ins_pattern"], ins_dict["up_seq"], local_aligner)
            dw_similar_score = cal_align_score(ins_dict["ins_pattern"], ins_dict["dw_seq"], local_aligner)
            if up_similar_score > 0.8 or dw_similar_score > 0.8:
                 ins_ls.append(ins_id)
        pd.DataFrame({"repeat_pattern_ins": ins_ls}).to_csv(output.ins_ls, index=False)
        
 
rule summary:
    input:
        rc_s=config["file"]["repeat_check"]["ins_ls"],
        ne_s=config["file"]["needle_summary"],
        replace_s=config["file"]["replace_repeatmasker"]
    output:
        config["file"]["decision_summary"]
    run:
        needle_df = pd.read_csv(config["file"]["needle_summary"], sep="\t")
        up_dup_df = needle_df[(needle_df["up_align_q_cov"]>=0.95) & (needle_df["up_align_q_cov"]<=1.05) \
                              & (needle_df["up_iden_perc"]>=0.9) & (needle_df["up_end_dis"]<=5)]
        dw_dup_df = needle_df[(needle_df["dw_align_q_cov"]>=0.95) & (needle_df["dw_align_q_cov"]<=1.05) \
                              & (needle_df["dw_iden_perc"]>=0.9) & (needle_df["dw_start_dis"]<=5)]
        dup_ls = set(up_dup_df["Unnamed: 0"]).union(set(dw_dup_df["Unnamed: 0"]))
        
        info_df = pd.read_csv(config["file"]["replace_repeatmasker"], sep="\t")
        info_df.set_index("Unnamed: 0", inplace=True)
        info_df["decision_type"] = info_df["final_decision"].apply(lambda x: x=="Unmask" and "Unmask" or ",".join(set(["_".join(i.split("_")[2:]).split("/")[0] for i in x.strip().split("|")])))
        
        pattern_ls = ["SINE","LINE","LTR","DNA","Retroposon","RNA","Simple_repeat","Tandem_repeat", "Satellite", "Low_complexity"]
        te_ls = ["SINE","LINE","LTR","DNA","Retroposon","RNA"]
        lowcomp_ls = ["Simple_repeat","Tandem_repeat", "Satellite"]
        
        lowcomp_idx_ls = []
        for idx, row in info_df.iterrows():
            if row["decision_type"] in lowcomp_ls:
                anno_len = sum([int(i.split("_")[1])-int(i.split("_")[0]) for i in row["final_decision"].split("|")])
                ## by anno ratio filtering
                if anno_len/row["total_len"] >= 0.9:
                    lowcomp_idx_ls.append(idx)
                ## by anno len filtering for short sequence
                if row["final_decision"].count("|")==0:
                    start, end = int(row["final_decision"].split("_")[0]), int(row["final_decision"].split("_")[1])
                    if (start <= 5 and row["total_len"]-end <=50) or (row["total_len"]-end <=5 and start <= 50):
                        lowcomp_idx_ls.append(idx)
        
        te_idx_ls = []
        for idx, row in info_df.iterrows():
            for s_value in row["decision_type"].split(","):
                if s_value in te_ls:
                    te_idx_ls.append(idx)

        pattern_df = pd.read_csv(config["file"]["repeat_check"]["ins_ls"])
        segment_dup = set(dup_ls).difference(set(pattern_df['repeat_pattern_ins'])).difference(set(lowcomp_idx_ls))
        repeat_dup = set(pattern_df['repeat_pattern_ins']).difference(set(te_idx_ls)).union(set(lowcomp_idx_ls))
        info_df.reset_index(inplace=True)
        info_df["segment_dup"] = info_df["Unnamed: 0"].apply(lambda x: x in segment_dup and 1 or 0)
        info_df["repeat_dup"] = info_df["Unnamed: 0"].apply(lambda x: x in repeat_dup and 1 or 0)
        info_df.to_csv(output[0], sep="\t", index=False)
        
