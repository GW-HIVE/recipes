"""
###
Full list of actions:
- union
- addconstantfield
- makecombofield
- filterin
- filterout
- expandrows
- transposecols
- addnormalizedcol
- addaveragecol
- splitcol
- jointables

###
union 
###
###Summary

Merge multiple tables and remove duplicate lines.

###Format:

{
	"id":#,
	"type":"output",
	"inputtables":[#,#,#],
	"action":{
		"name":"union"
	}
}

###Notes:
This action can take multiple inputtable tables. 

See use-case recipe example_part1.json for usage example.

###
addconstantfield 
###
###Summary

Add a new field to the table that has the same constant entry for each line.

###Format:

{
	"id":#,
	"type":"output",
	"inputtables":[#,#,#],
	"action":{
		"name":"addconstantfield",
		"newfields":{
			"example_field_1":"example_entry_1",
			"example_field_2":"example_entry_2"
		}
	}
}

###Notes: 
This action can take multiple inputtable tables. 

See use-case recipe example_part2.json for usage example.

###
makecombofield
###
###Summary

Combine two fields into one field. 

###Format:

{
	"id":#,
	"type":"output",
	"inputtables":[#],
	"action":{
		"name":"makecombofield",
		"fieldlist":[
			"example_field_1",
			"example_field_2"
		],
		"merge_char":"|",
		"combofield":"example_combo_field"
	}
}

###Notes: 
This action can take a single inputtable table. 

The above example makes the combined field example_combo_field from example_field_1 and example_field_2. 

The entries in each of the fields will merge into the combined field and be separated by the merge character, in the example the charcater is |.


See use-case recipe example_part1.json for usage example.

###
filterin
###
###Summary

Filter only certain lines from one table into another.

Needs review from Robel.

###Format:

{
	"id":#,
	"type":"output",
	"inputtables":[#],
	"action":{
		"name":"filterin",
		"operation":"AND",
		"conditionlist":[
			"field":"example_field_1",
			"value":[
				"ABC",
				"-",
				""
			],
			"operation":"in"
	}
}

###Notes: 
This action can take a single inputtable table. 

In the above example, only the lines with the values ABC, -, or an empty field will be included in to a new table.

See use-case recipe example_part2.json for usage example.

Question for Robel:
What is the purpose of "operation"?

###
filterout
###
###Summary

Filter out only certain lines from a table to create a table without those entries.

Needs review from Robel.

###Format:

{
	"id":#,
	"type":"output",
	"inputtables":[#],
	"action":{
		"name":"filterout",
		"operation":"AND",
		"conditionlist":[
			"field":"example_field_1",
			"value":[
				"DEF"
				"123"
				"456"
			],
			"operation":"in"
	}
}

###Notes: 
This action can take a single inputtable table. 

In the above example, lines with the values DEF, 123, or 456 in the field example_field_1 will be excluded from a new table, and all other lines with other values for that field will be included.

See use-case recipe example_part1.json for usage example.

Question for Robel:
What is the purpose of "operation"?

###
expandrows
###
###Summary

For a single line that has a field with multiple data points in one field, create multiple lines where each line has a single data point from that field. 

###Format:

{
	"id":#,
	"type":"output",
	"inputtables":[#],
	"action":{
		"name":"expandrows",
		"expansionfield":"example_field_1",
		"expansiondelim":";"
		}
	}
}

###Notes: 
This action can take a single inputtable table. 

In the above example, the field exmple_field_1 will be used to expand rows. For every line that has the expansion delim ;, the data in that row will be separated and a new line created for each data point.

An example:

Table input: 
"example_field_1"
"data_point_1;data_point_2","information_abc"

Table output from expandrows above: 
"example_field_1","example_field_2"
"data_point_1","information_abc"
"data_point_2","information_abc"

See use-case recipe example_part2.json for usage example.

###
transposecols
###
###Summary

The rows switch with columns and vice versa(?)

Needs review from Robel.

###Format:

{
	"id":#,
	"type":"output",
	"inputtables":[#],
	"action":{
		"name":"transposecols",
		"startcolidx":1,
		"newfieldone":"new_example_field_1"
		"newfieldtwo":"new_example_field_2"
		}
	}
}

###Notes: 
This action can take a single inputtable table.

See use-case recipe example_part2.json for usage example.

Question for Robel:
How is this action used?

###
addnormalizedcol
###
###Summary

Create a new column that takes the numerical data from an existing column and normalizes the numbers according to (?)

Needs review from Robel.

###Format:

{
	"id":#,
	"type":"output",
	"inputtables":[#],
	"action":{
		"name":"addnormalizedcol",
		"startcolidx":1,
		"newfield":"normalized_example_field_1"
		"anchorfields":["example_field_2"]
		}
	}
}

###Notes: 
This action can take a single inputtable table.

See use-case recipe example_part2.json for usage example.

Question for Robel:
How is normalization calculated?

###
addaveragecol
###
###Summary

Create a new column that takes the numerical data from an existing column and averages the numbers according to (?)

Needs review from Robel.

###Format:

{
	"id":#,
	"type":"output",
	"inputtables":[#],
	"action":{
		"name":"addaveragecol",
		"field":"example_field_1",
		"newfield":"average_example_field_1"
		"anchorfields":["example_field_2"]
		}
	}
}

###Notes: 
This action can take a single inputtable table.

See use-case recipe example_part2.json for usage example.

Question for Robel:
How is the average being calculated?


###
splitcol
###
###Summary

Split a single column into two separate columns after specifying a deliminating character that separates the values. 

###Format:

{
	"id":#,
	"type":"output",
	"inputtables":[#],
	"action":{
		"name":"splitcol",
		"field":"example_field_1",
		"newfields":["example_field_1a","example_field_1b"],
		"delim":"|"
		}
	}
}

###Notes: 
This action can take a single inputtable table.

See use-case recipe example_part1.json for usage example.


###
jointables
###
###Summary

Create a new table that uses an anchor field to combine columns from multiple tables. Data with the same anchor field value will be added to the same row in the new field. 

###Format:

{
	"id":#,
	"type":"output",
	"inputtables":[#,#,#],
	"action":{
		"name":"jointables",
		"anchorfields":["example_field_1","example_field_2"]
		}
	}
}

###Notes: 
This action can take multiple inputtable table.

See use-case recipe example_part1.json for usage example.
"""

import os,sys
import json
import csv
#import requests

from optparse import OptionParser
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist
from Bio.pairwise2 import format_alignment

from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

import commands
import glob
import datetime
import pytz

sys.path.append('../lib/')
import csvutil


__version__="1.0"
__status__ = "Dev"

def dump_log(log_file, msg, flag):

    ts = datetime.datetime.now(pytz.timezone('US/Eastern')).strftime('%Y-%m-%d %H:%M:%S %Z%z')
    if flag == "write":
        with open(log_file, "w") as FW:
            FW.write("%s\t%s\n" % (ts,msg))
    else:
        with open(log_file, "a") as FA:
            FA.write("%s\t%s\n" % (ts, msg))

    return


def execute_recipe(dataset):

    masterlist_file = "unreviewed/biomarker_masterlist.csv"
    recipe_file = "recipes/%s.json" % (dataset)
    log_file = "logs/%s.log" % (dataset)
    
    recipe_json = json.loads(open(recipe_file, "r").read())
    
    #make sure one of the output files is named as dataset name
    flag = False
    for obj in recipe_json["outputfiles"]:
        if obj["filepath"].find(dataset) != -1:
            flag = True
    if flag == False:
        print "None of the output file names is the same as the dataset (%s)" % (dataset)
        sys.exit()



    dump_log(log_file, "started making dataset", "write")

    df_list = []
    for obj in recipe_json["tablelist"]:
        if obj["type"] != "input":
            continue
        data_frame = {}
        sep = obj["separator"]
        file_type = obj["filepath"].split(".")[-1]
        if file_type in ["csv", "tsv"]:
            if "sizecategory" in obj:
                if obj["sizecategory"] in ["regular"]:
                    csvutil.load_sheet(data_frame, obj["filepath"], obj["fields"], sep)
                else:
                    csvutil.load_large_sheet(data_frame, obj["filepath"], obj["fields"], sep)
            else:
                csvutil.load_large_sheet(data_frame, obj["filepath"], obj["fields"], sep)
        elif file_type in ["json"]:
            csvutil.load_sheet_from_json(data_frame, obj["filepath"], obj["fields"])
        df_list.append(data_frame)

    dump_log(log_file, "finished loading input", "append")


    for obj in recipe_json["tablelist"]:
        if obj["type"] != "output":
            continue
        if "action" in obj:
            o = obj["action"]
            if o["name"] == "union":
                data_frame = csvutil.union_tables(df_list,obj["inputtables"])
                df_list.append(data_frame)
            elif o["name"] == "addconstantfields":
                table_ind = obj["inputtables"][0] - 1
                data_frame = csvutil.add_constant_fields(df_list,table_ind, o["newfields"])
                df_list.append(data_frame)
            elif o["name"] == "makecombofield":
                table_ind = obj["inputtables"][0] - 1
                data_frame = csvutil.add_combo_field(df_list,table_ind,o)
                df_list.append(data_frame)
            elif o["name"] == "filterin":
                table_ind = obj["inputtables"][0] - 1
                data_frame = csvutil.filter_in_records(df_list,table_ind,o)
                df_list.append(data_frame)
            elif o["name"] == "filterout":
                table_ind = obj["inputtables"][0] - 1
                data_frame = csvutil.filter_out_records(df_list,table_ind,o)
                df_list.append(data_frame)
            elif o["name"] == "expandrows":
                table_ind = obj["inputtables"][0] - 1
                data_frame = csvutil.expand_rows(df_list,table_ind,o)
                df_list.append(data_frame)
            elif o["name"] == "transposecols":
                table_ind = obj["inputtables"][0] - 1
                data_frame = csvutil.transpose_cols(df_list,table_ind,o)
                df_list.append(data_frame)
            elif o["name"] == "addnormalizedcol":
                table_ind = obj["inputtables"][0] - 1
                data_frame = csvutil.add_normalized_col(df_list,table_ind,o)
                df_list.append(data_frame)
            elif o["name"] == "addaveragecol":
                table_ind = obj["inputtables"][0] - 1
                data_frame = csvutil.add_average_col(df_list,table_ind,o)
                df_list.append(data_frame)
            elif o["name"] == "splitcol":
                table_ind = obj["inputtables"][0] - 1
                data_frame = csvutil.split_col(df_list,table_ind,o)
                df_list.append(data_frame)
            elif o["name"] == "jointables":
                table_ind_one = obj["inputtables"][0] - 1
                table_ind_two = obj["inputtables"][1] - 1
                data_frame = csvutil.join_tables(df_list[table_ind_one], df_list[table_ind_two],
                        o["anchorfields"])
                df_list.append(data_frame)
             



    dump_log(log_file, "finished performing actions", "append")


    for obj in recipe_json["outputfiles"]:
        #populate id master list if masterlistfilterfield is set
        masterlist_dict = {}
        masterlist_field = ""
        if "masterlistfilterfield" in obj:
            if obj["masterlistfilterfield"] != "":
                masterlist_field = obj["masterlistfilterfield"]
                data_frame = {}
                sep = "\",\""
                fields = [masterlist_field]
                csvutil.load_large_sheet(data_frame, masterlist_file, fields, sep)
                for row in data_frame["data"]:
                    masterlist_dict[row[0].upper()] = True
        
        table_ind = obj["tableid"] - 1
        out_file = obj["filepath"]
        data_frame = df_list[table_ind]

        FW = open(out_file, "w")
        newrow = []
        field_map = {}
        f_list_new = []
        for o in obj["fields"]:
            field_map[o["name"]] = o["newname"]
            field_map[o["newname"]] = o["name"]
            f_list_new.append(o["newname"])
        FW.write("\"%s\"\n" % ("\",\"".join(f_list_new)))
       
        dump_log(log_file, "started processing output rows", "append")
        row_count = 0
        for row in data_frame["data"]:
            row_count += 1
            if row_count%100000 == 0:
                dump_log(log_file, "processed %s output rows" % (row_count), "append")
            
            newrow = []
            for f_new in f_list_new:
                f_old = field_map[f_new]
                value = row[data_frame["fields"].index(f_old)]
                if f_new in ["log2fc"]:
                    value = "%1.2f" %(float(value))
                elif f_new in ["pvalue","adjpvalue"]:
                    value = "%1.2e" %(float(value))
                newrow.append(value)
            
            #filter out if masterlist_field is set
            if masterlist_field != "":
                masterlist_field_old = field_map[masterlist_field]
                masterlist_value = row[data_frame["fields"].index(masterlist_field_old)].upper()
                #print masterlist_field_old, masterlist_field, masterlist_value, masterlist_value in masterlist_dict 
                if masterlist_value not in masterlist_dict:
                    
                    continue
            FW.write("\"%s\"\n" % ("\",\"".join(newrow)))
        FW.close()

    dump_log(log_file, "finished making dataset", "append")

    return





###############################
def main():

    usage = "\n%prog  [options]"
    parser = OptionParser(usage,version="%prog " + __version__)
    parser.add_option("-d","--dataset",action="store",dest="dataset",help="[idmapping, transcriptlocus]")


    (options,args) = parser.parse_args()
    for file in ([options.dataset]):
        if not (file):
            parser.print_help()
            sys.exit(0)

    global config_obj
    global species_obj

    dataset = options.dataset

    config_obj = json.loads(open("../conf/config.json", "r").read())
    species_obj = config_obj["speciesinfo"]
 
    global path_obj
    path_obj = config_obj["pathinfo"]

    
    execute_recipe(dataset)
    return






if __name__ == '__main__':
        main()

