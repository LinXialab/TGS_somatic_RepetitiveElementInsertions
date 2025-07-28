rule summary_trf_sdust:
    input:
        config["file"]["trf_size_check"]
    output:
        config["file"]["sdust_trf_summary"]
    run:
        result_dict = dict()
        for ins_id in sv_ls:
            sdust_path = "{}/{}.out".format(config["path"]["output_sdust"], ins_id)
            trf_path = "{}/{}.out".format(config["path"]["output_trf"], ins_id)
            fasta_path = get_fasta_path(ins_id)
            records = SeqIO.parse(fasta_path, "fasta")
            record = list(records)[0]
            total_len = len(record.seq)
            trf_fasta_path = "{}/{}.fasta".format(config["path"]["pattern_f"], ins_id)
            if not os.path.isfile(trf_fasta_path):
                trf_seq, trf_len, trf_cor = None, 0, []
            else:
                trf_records = SeqIO.parse(trf_fasta_path, "fasta")
                trf_seq = ",".join([str(i.seq) for i in list(trf_records)])
                trf_len, trf_cor = trf_parser(trf_path)
            if os.path.getsize(sdust_path) == 0:
                sdust_len = 0
                sdust_cor = []
            else:
                df = pd.read_csv(sdust_path, sep="\t", names=["id", "start", "end"])
                df["len"] = df["end"] - df["start"]
                df["cor"] = df.apply(lambda x: "{}_{}".format(x.start, x.end), axis=1)
                sdust_len = df["len"].sum()
                sdust_cor = "|".join(df["cor"].tolist())
            result_dict.setdefault("sdust_ratio", {}).setdefault(ins_id, sdust_len/total_len)
            result_dict.setdefault("sdust_cor", {}).setdefault(ins_id, sdust_cor)
            result_dict.setdefault("trf_ratio", {}).setdefault(ins_id, "{:.4f}".format(trf_len/total_len))
            result_dict.setdefault("trf_cor", {}).setdefault(ins_id, "|".join(["{}_{}".format(i[0], i[1]) for i in trf_cor]))
            result_dict.setdefault("trf_seq", {}).setdefault(ins_id, trf_seq)
            result_dict.setdefault("total_len", {}).setdefault(ins_id, total_len)
        pd.DataFrame(result_dict).to_csv(output[0], sep="\t")

rule replace_repeatmasker:
    input:
        tf_s=rules.summary_trf_sdust.output,
        rm_s=rules.summary_repeatmasker.output,
    output:
        config["file"]["replace_repeatmasker"] 
    run:
        lowcomp_df = pd.read_csv(config["file"]["sdust_trf_summary"], sep="\t")
        lowcomp_df.set_index("Unnamed: 0", inplace=True)
        repeat_df = pd.read_csv(config["file"]["repeatmasker_summary"], sep="\t")
        repeat_df.set_index("Unnamed: 0", inplace=True)
        miss_ls = insert_df.index.difference(lowcomp_df.index)
        
        info_df = lowcomp_df.join(repeat_df)
        info_df.fillna(0, inplace=True)
        info_df["repeatMasker_type"] = info_df["len"].apply(lambda x: x in ["0",0] and "Unmask" or ",".join(set(["_".join(i.split("_")[1:]) for i in x.strip().split("|")])))
        info_df["repeatMasker_cat"] = info_df["repeatMasker_type"].apply(get_te_cat)
        info_df["sdust_cor"] = info_df.apply(concat_sdust_type, axis=1)
        info_df["trf_cor"] = info_df.apply(lambda x: concat_trf_type(x), axis=1)
        
        ## make decision
        decision = []
        replaced = []
        for idx, row in info_df.iterrows():
            cat = row["repeatMasker_cat"]
            repeat_type = row["repeatMasker_type"]
            ## add trf anno only if repeatmasker TE-unannoted length larger than trf annotation length
            if "TE" in cat:
                type_s, change = cross_check_cor_te(row["coordinate"], row["trf_cor"])
                decision.append(type_s)
                replaced.append(change)
            elif cat == "Unmask":
                decision.append(row["trf_cor"] and row["trf_cor"] or (row["sdust_cor"] and row["sdust_cor"] or "Unmask"))
                replaced.append(row["trf_cor"] and 1 or (row["sdust_cor"] and 1 or 0))
            ## simple_repeat, low_complexity, satellite
            else:
                type_s, change = cross_check_cor_te(row["coordinate"], row["trf_cor"])
                decision.append(type_s)
                replaced.append(change)
        
        info_df["final_decision"] = decision
        info_df["replaced"] = replaced
        info_df.to_csv(output[0], sep="\t")
