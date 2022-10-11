rule get_refs:
    output:
        expand(
            "references/{ref}",
            ref=config["references"].values()
        )
    run:
        print(config["references"])
        for ref in config["references"].values():
            print(ref)
            shell("gsutil -o \"GSUtil:parallel_process_count=1\" -m cp gs://genomics-public-data/resources/broad/hg38/v0/{ref} references/", ref=ref)
