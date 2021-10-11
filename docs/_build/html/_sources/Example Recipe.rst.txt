Example Recipe
---------------

*Here is the example.json recipe file with comments*::

    {
        "outputfiles": [
            {
                "fields": [
                
                    # for fields, the name must be a name that has been given in the tablelist section below
                    # these are often field names from the input datasets, but also can be fieldnames generated from actions
                    # newname is the fieldname that will be generated in the output file and replace the old field name
                    {
                        "newname": "gene_symbol",
                        "name": "Gene Symbol"
                    },
            {
                        "newname": "uniprot_ac",
                        "name": "Uniprot AC"
                    },
            {
                        "newname": "chr_num",
                        "name": "Chr Num"
                    },
            {
                        "newname": "chr_pos",
                        "name": "Chr Pos"
                    },
            {
                        "newname": "ref_nt",
                        "name": "Ref NT"
                    },
            {
                        "newname": "alt_nt",
                        "name": "Alt NT"
                    }
                ], 
                
            # the tableid denotes what table from the tablelist will be used for this output file
                "tableid": 3, 
                
            # the masterlistfilter field denotes the fieldname that is used as the primary field name for filtering
                "masterlistfilterfield":"gene_symbol",
                
            # the filepath gives the name and path for the outputfile
                "filepath": "unreviewed/example_dataset.csv"
            }
        ], 
        "tablelist": [
            {
            # every table in tablelist has a tableid, which can be referenced by other acgtions or the outputfile 
                "id": 1,
                
            # input denotes we are taking a file in from the filepath listed
            "type": "input",
                
            # the separator denotes the structure of the input file
            # here the separator is 'doublequote comma doublequote' 
            # so each line of the input file looks like this "a","b","c" 
            "separator": "\",\"",
                
            # list of the fields in the inout file
                "fields": [
                    "Uniprot AC",
                    "Chr Num", 
                    "Chr Pos",
                    "Ref NT", 
                    "Alt NT",
            "EXAC Alt Frequency",
            "1000 Genomes Alt Frequency"
                ],
                "filepath": "downloads/example_mutations.csv"
            },
        {
            # same as above for table '"id":1', but for a different input file
            "id": 2,
            "type": "input",
                "separator": "\",\"", 
                "fields": [
                    "uniprot accession",
                    "Gene Symbol" 
                ],  
                "filepath": "downloads/example_mapping.csv"
            },
        {
            # this table performs an action
            "id": 3,
                
            # the type is output so it will be generating a new table based on the action used
            "type": "output",  
                
            # the tables that are input to generate this new table. 
            # in this case our input files from tables 1 and 2
            "inputtables":[1,2],
                
            # what action will be taken to make the new table
            "action":{
                # the name gives the specific action to take
            "name":"jointables",
                    
            # different actions have their own parameters. 
            # in this case the jointables action has only the parameter "anchorfields".
            # specific details on actions are below in the README
            "anchorfields":["Uniprot AC", "uniprot accession"]
                },
            }
        ]
    }
