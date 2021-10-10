"""

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

