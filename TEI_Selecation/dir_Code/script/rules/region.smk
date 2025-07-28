rule determine_region:
    input:
        summary_dict=config["file"]["repeat_check"]["summary_dict"],
        cor_dict=config["file"]["repeat_check"]["cor_dict"],
        ne_s=config["file"]["needle_summary"],
        decision_s=config["file"]["decision_summary"]
    output:
        config["file"]["final_summary"]
    run:
        info_df = pd.read_csv(input.decision_s, sep="\t")
        info_df.set_index("Unnamed: 0", inplace=True)
        needle_df = pd.read_csv(input.ne_s, sep="\t")
        needle_df.set_index("Unnamed: 0", inplace=True)
        segment_dup_ls = info_df[info_df["segment_dup"] == 1].index.tolist()
       
        ## set default dup region
        dup_region_d_ls = []
        for idx, row in info_df.iterrows():
            ins_region = insert_df.loc[idx, "ins_loc_hg38"]
            ref_chr, ref_start = ins_region.split("_")[0], int(ins_region.split("_")[1])
            ref_end = ref_start+1
            dup_region_d_ls.append("{}_{}_{}".format(ref_chr, ref_start-10, ref_end+10))
        info_df["dup_region"] = dup_region_d_ls

        ## determine segment dup dup region
        dup_region_dict = dict()
        for idx, row in info_df.loc[segment_dup_ls].iterrows():
            ins_region = insert_df.loc[idx, "ins_loc_hg38"]
            ref_chr, ref_start = ins_region.split("_")[0], int(ins_region.split("_")[1])
            ref_end = ref_start+1
            ins_len = int(row["total_len"])-100
            seg_info = needle_df.loc[idx]
            up_qs_dis = seg_info["up_start_dis"]
            dw_qs_dis = seg_info["dw_end_dis"]
            if abs(seg_info["up_end_dis"]) <= 5 and seg_info["up_align_q_cov"]>=0.95 \
               and seg_info["up_align_q_cov"]<=1.05 and seg_info["up_iden_perc"]>=0.9:
                up_s_start = ref_start+1-(ins_len+50)
                region_start = up_qs_dis > 0 and up_s_start+up_qs_dis or up_s_start
                region_end = ref_start
            else:
                dw_s_end = ref_end+ins_len+50
                region_start = ref_end
                region_end = dw_qs_dis <0 and dw_s_end+dw_qs_dis or dw_s_end
            dup_region = "{}_{:.0f}_{:.0f}".format(ref_chr, region_start, region_end)
            dup_region_dict[idx] = dup_region
            
        ## determine repeat dup dup region by pattern pairwise alignemnt
        repeat_dup_ls = info_df[info_df["repeat_dup"] == 1].index.tolist()
        with open(input.cor_dict, "rb") as handle:
            ins_cor_dict = pickle.load(handle)
        with open(input.summary_dict, "rb") as handle:
            ins_pattern_dict = pickle.load(handle)
        ### check repeat dup with pattern alignment information
        repeat_cor_ls = set(repeat_dup_ls).intersection(ins_cor_dict.keys())
        repeat_info_dict = dict()
        for ins_id in repeat_cor_ls:
            ins_region = insert_df.loc[ins_id, "ins_loc_hg38"]
            ref_chr, ref_start = ins_region.split("_")[0], int(ins_region.split("_")[1])
            ref_end = ref_start+1
            if "up_s_align_cor" in ins_cor_dict[ins_id]:
                up_seq_len = len(ins_pattern_dict[ins_id]["up_seq"])*ins_cor_dict[ins_id]["up_s_times"]
                up_align_start = ins_cor_dict[ins_id]["up_s_align_cor"][0]
                up_ref_start = ref_start - (up_seq_len-up_align_start)
                pattern = "({}){}".format(ins_pattern_dict[ins_id]["ins_pattern"], ins_cor_dict[ins_id]["up_q_times"])
                repeat_region = "{}_{}_{}".format(ref_chr, up_ref_start, ref_start)
                repeat_info_dict.setdefault(ins_id, {}).update({"region": repeat_region, "pattern": pattern})
            else:
                dw_align_end = ins_cor_dict[ins_id]["dw_s_align_cor"][1]
                dw_ref_end = ref_end + dw_align_end
                pattern = "({}){}".format(ins_pattern_dict[ins_id]["ins_pattern"], ins_cor_dict[ins_id]["dw_q_times"])
                repeat_region = "{}_{}_{}".format(ref_chr, ref_end, dw_ref_end)
                repeat_info_dict.setdefault(ins_id, {}).update({"region": repeat_region, "pattern": pattern})
        repeat_region_dict = dict()
        for idx, aa_dict in repeat_info_dict.items():
            repeat_region_dict[idx]=aa_dict["region"]
        
        ### check repeat dup without pattern alignment information, e.g. low complexity
        anno_ls = set(repeat_dup_ls).difference(ins_cor_dict.keys())
        anno_repeat_region_dict = dict()
        for idx, row in info_df.loc[anno_ls].iterrows():
            ins_region = insert_df.loc[idx, "ins_loc_hg38"]
            ref_chr, ref_start = ins_region.split("_")[0], int(ins_region.split("_")[1])
            ref_end = ref_start+1
            ins_len = row["total_len"]
            if row["trf_seq"] != "0" and row["trf_seq"].count(",")==0:
                pattern_len = len(row["trf_seq"])
                include_cor_ls = []
                for cor in row["trf_cor"].split("|"):
                    start, end = map(int, cor.split("_")[:2])
                    if end <= 50:
                        continue
                    elif ins_len - start <= 50:
                        continue
                    else:
                        include_cor_ls.append((start, end))
                if not include_cor_ls:
                    continue
                start, end = include_cor_ls[0][0], include_cor_ls[0][1]
                repeat_start, repeat_end, ref_repeat_times = get_repeat_cor(pattern_len, start, end, ins_len, ref_start, ref_end, ref_chr)
            else:
                if "Satellite" in row["final_decision"] or row["trf_seq"]=="0":
                    include_cor = row["final_decision"].split("_")
                    start, end = int(include_cor[0]), int(include_cor[1])
                    ref_repeat_times = 1
                    repeat_start = ref_start - (51-start)
                    repeat_end = ref_start + 49-(ins_len-end)
                else:
                    anno_cor_ls = row["trf_cor"].split("|")
                    if int(anno_cor_ls[0].split("_")[0]) > 10:
                        start, end = map(int, anno_cor_ls[-1].split("_")[:2])
                        pattern_len = len(row["trf_seq"].split(",")[-1])
                        repeat_start, repeat_end, ref_repeat_times = get_repeat_cor(pattern_len, start, end, ins_len, ref_start, ref_end, ref_chr)
                    else:
                        start, end = map(int, anno_cor_ls[0].split("_")[:2])
                        pattern_len = len(row["trf_seq"].split(",")[0])
                        repeat_start, repeat_end, ref_repeat_times = get_repeat_cor(pattern_len, start, end, ins_len, ref_start, ref_end, ref_chr)
            repeat_region = "{}_{}_{}".format(ref_chr, repeat_start, repeat_end)
            anno_repeat_region_dict[idx] = repeat_region

        dup_region_dict.update(repeat_region_dict)
        dup_region_dict.update(anno_repeat_region_dict)
        info_df.reset_index(inplace=True)
        info_df["dup_region"] = info_df.apply(lambda x: x["Unnamed: 0"] in dup_region_dict.keys() and dup_region_dict[x["Unnamed: 0"]] or x["dup_region"], axis=1)
        info_df["clean_decision"] = info_df.apply(clean_anno, axis=1)
        info_df.rename(columns={"Unnamed: 0": "ins_id"}, inplace=True)
        info_df.to_csv(output[0], sep="\t", index=False)

