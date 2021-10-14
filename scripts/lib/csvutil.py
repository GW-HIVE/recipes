'''
This script provides performs given actions on datasets, as specified by a recipe json.

Full list of documented actions:
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

Input:
#######
    * df_list : The data tables taken as input in the recipe json.
    * table_ind : The data tables created as output in the recipe json.
    * action_obj : The specific action performed in the recipe json.

Output:
#######
    * data_frame : The final data table produced by a given action.

'''

import csv
import json


def split_col(df_list, table_ind, action_obj):
'''
###
Split Column
###
###Summary

Split a single column into two separate columns.

Input:
#######
    * field : The given field to split.
    * newfields : The name of the two new fields.
    * delim : The deliminator character where the value will be split.

###Recipe Action Format:

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
'''
    field_list = df_list[table_ind]["fields"]
    data_frame = {"fields":field_list, "data":df_list[table_ind]["data"]}
   
    field_idx = field_list.index(action_obj["field"])
    delim = action_obj["delim"]
    new_field_list = field_list[0:field_idx] + action_obj["newfields"] 
    new_field_list += field_list[field_idx+1:]

    row_list = []
    for row in data_frame["data"]:
        val_list = row[field_idx]
        newrow = row[0:field_idx] + row[field_idx].split(delim) + row[field_idx+1:]
        row_list.append(newrow)
    
    data_frame["data"] = row_list
    data_frame["fields"] = new_field_list

    return data_frame



def add_average_col(df_list, table_ind, action_obj):
'''
###
Add Average Column
###
###Summary

Create a new column that takes the numerical data from an existing column and averages the numbers according to (?)

####Average calculation

x_sum/x_count

x_sum = sum of all values for given field
x_count = number of values for a given field

Input:
#######
    * field : The given field to average.
    * newfield : The name of the field to display the average.
    * anchorfields : A list of the fields that represent unique rows.

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
'''

    field_list = df_list[table_ind]["fields"]
    data_frame = {"fields":field_list, "data":df_list[table_ind]["data"]}
    
    row_grps = {}
    for row in data_frame["data"]:
        key_list = []
        for f in action_obj["anchorfields"]:
            key_list.append(row[field_list.index(f)])
        group_key = " ".join(key_list)
        if group_key not in row_grps:
            row_grps[group_key] = []
        row_grps[group_key].append(row)

    row_list = []
    for group_key in row_grps:
        x_count, x_sum, x_max,x_min = 0.0, 0.0, -100000.0, 1000000.0
        for row in row_grps[group_key]:
            x = float(row[field_list.index(action_obj["field"])])
            x_sum += x
            x_count += 1.0
        x_average = "%.5f" % (x_sum/x_count)
        for row in row_grps[group_key]:
            newrow =  row + [x_average]
            row_list.append(newrow)
    
    data_frame["data"] = row_list
    data_frame["fields"] = field_list + [action_obj["newfield"]]
    return data_frame



    

def add_normalized_col(df_list, table_ind, action_obj):
'''
###
Add Normalized Column
###
###Summary

Create a new column with normalized values from a given column.

####Normalization calculation

x/x_sum

x = value in given field and per row
x_sum = sum of all values for given field

Input:
#######
    * field : The given field to normalize
    * newfield : The name of the new normalized field
    * anchorfields : A list of the fields that represent unique rows

###Format:

{
	"id":#,
	"type":"output",
	"inputtables":[#],
	"action":{
		"name":"addnormalizedcol",
		"field":"example_field_1",
		"newfield":"normalized_example_field_1"
		"anchorfields":["example_field_2"]
		}
	}
}
'''

    field_list = df_list[table_ind]["fields"]
    data_frame = {"fields":field_list, "data":df_list[table_ind]["data"]}
    
    row_grps = {}
    for row in data_frame["data"]:
        key_list = []
        for f in action_obj["anchorfields"]:
            key_list.append(row[field_list.index(f)])
        group_key = " ".join(key_list)
        if group_key not in row_grps:
            row_grps[group_key] = []
        row_grps[group_key].append(row)

    row_list = []
    for group_key in row_grps:
        x_sum, x_max,x_min = 0.0, -100000.0, 1000000.0
        for row in row_grps[group_key]:
            x = float(row[field_list.index(action_obj["field"])])
            x_max = x if x > x_max else x_max
            x_min = x if x < x_min else x_min
            x_sum += x

        for row in row_grps[group_key]:
            x = float(row[field_list.index(action_obj["field"])])
            #x_normalized = (x-x_min)/(x_max-x_min)
            x_normalized = "%.5f" % (x/x_sum)
            newrow =  row + [x_normalized]
            row_list.append(newrow)


    data_frame["data"] = row_list
    data_frame["fields"] = field_list + [action_obj["newfield"]]
    return data_frame



def transpose_cols(df_list, table_ind, action_obj):
'''
###
Transpose Columns
###
###Summary

Rows become columns and columns become rows.

Input:
#######
    * startcolidx : The left most column to begin transposing into rows.
    * newfieldone : A new field to add after the transposition.
    * newfieldtwo : A second new field to add after the transposition.

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
'''


    field_list = df_list[table_ind]["fields"]
    data_frame = {"fields":field_list, "data":df_list[table_ind]["data"]}
    start_col_idx = action_obj["startcolidx"]
    new_field_one = action_obj["newfieldone"]
    new_field_two = action_obj["newfieldtwo"]

    row_list = []
    field_list_new = field_list[0:start_col_idx] + [new_field_one, new_field_two]
    for row in data_frame["data"]:
        for j in xrange(start_col_idx, len(row)):
            v = row[j]
            newrow = row[0:start_col_idx] + [field_list[j],v]
            row_list.append(newrow)

    data_frame["data"] = row_list
    data_frame["fields"] = field_list_new

    return data_frame


def expand_rows(df_list, table_ind, action_obj):
'''
###
Expand Rows
###
###Summary

Create new rows from a single field value that contains multiple entries.

Input:
#######
    * expansionfield : The field containing multiple entries to expand into new rows.
    * expansiondelim: The deliminating character that separates the entries.

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

'''

    field_list = df_list[table_ind]["fields"]
    data_frame = {"fields":field_list, "data":df_list[table_ind]["data"]}
    expansion_field = action_obj["expansionfield"]
    expansion_delim = action_obj["expansiondelim"]
    row_list = []
    for row in data_frame["data"]:
        field_idx = field_list.index(expansion_field)
        val_list = list(set(row[field_idx].strip().split(expansion_delim)))
        
        for val in val_list:
            tmprow = []
            for j in xrange(0, len(row)):
                tmprow.append(row[j] if j != field_idx else val)
            row_list.append(tmprow)
    
    data_frame["data"] = row_list

    return data_frame


def filter_out_records(df_list, table_ind, action_obj):
'''
###
Filter Out
###
###Summary

Filter out certain rows from a table.

Input:
#######
    * conditionlist : A list of further inputs to create the filter.
    * field : The field to apply the filter to.
    * value : A list of values in the field to filter out.
    * operation : Must be 'in' for filter to function.

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
'''
    field_list = df_list[table_ind]["fields"]
    data_frame = {"fields":field_list, "data":df_list[table_ind]["data"]}

    row_list = []
    for row in data_frame["data"]:
        flag_list = []
        for obj in action_obj["conditionlist"]:
            cond_field = obj["field"]
            cond_value = obj["value"]
            data_value = row[field_list.index(cond_field)]
            if obj["operation"] == "in":
                flag_list.append(data_value in cond_value)

        flag_list_unique = list(set(flag_list))
        if len(flag_list) == len(action_obj["conditionlist"]) and flag_list_unique == [True]:
            continue
        row_list.append(row)

    data_frame["data"] = row_list
    return data_frame


def filter_in_records(df_list, table_ind, action_obj):
'''
###
Filter In
###
###Summary

Include only values specified in a filter to a new table.

Input:
#######
    * conditionlist : A list of further inputs to create the filter.
    * field : The field to apply the filter to.
    * value : A list of values in the field to filter in.
    * operation : Must be 'in' for filter to function.

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
'''

    field_list = df_list[table_ind]["fields"]
    data_frame = {"fields":field_list, "data":df_list[table_ind]["data"]}

    row_list = []
    for row in data_frame["data"]:
        flag_list = []
        for obj in action_obj["conditionlist"]:
            cond_field = obj["field"]
            cond_value = obj["value"]
            data_value = row[field_list.index(cond_field)]
            if obj["operation"] == "in":
                flag_list.append(data_value in cond_value)

        flag_list_unique = list(set(flag_list))
        if len(flag_list) == len(action_obj["conditionlist"]) and flag_list_unique == [True]:
            row_list.append(row)

    data_frame["data"] = row_list
    return data_frame



def add_combo_field(df_list, table_ind, action_obj):
'''
###
Make Combo Field
###
###Summary

Combine two fields into one field.

Input:
#######
    * fieldlist : A list of fields to combine.
    * merge_char : The deliminating characater that will separate values from the combined field.
    * combofield : The name of the new combined field.

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
'''
    field_list = df_list[table_ind]["fields"] + [action_obj["combofield"]]
    data_frame = {"fields":field_list, "data":df_list[table_ind]["data"]}
   

    for row in data_frame["data"]:
        extra_row = []
        for f in action_obj["fieldlist"]:
            extra_row.append(row[field_list.index(f)].replace("\"", "").strip())
        row += [action_obj["merge_char"].join(extra_row)]

    return data_frame


def add_constant_fields(df_list, table_ind, new_data):
'''
###
Add Constant Field
###
###Summary

Add new fields to a table that each have a constant value per row.

Input:
#######
    * newfields : A list of fields and values to add.

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
'''

    new_data_field_list = new_data.keys()
    field_list = df_list[table_ind]["fields"] + new_data_field_list
    data_frame = {"fields":field_list, "data":df_list[table_ind]["data"]}

    extra_row = []
    for f in new_data_field_list:
        extra_row.append(new_data[f])
    for row in data_frame["data"]:
        row += extra_row

    return data_frame





def union_tables(df_list,ind_list):
'''
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
'''
    out_tbl = []
    first_ind = ind_list[0] - 1
    for i in ind_list:
        ind = i - 1
        out_tbl += df_list[ind]["data"]
    
    seen = {}
    data_frame = {"fields":df_list[first_ind]["fields"], "data":[]}
    for row in out_tbl:
        row_str = json.dumps(row)
        if row_str in seen:
            continue
        data_frame["data"].append(row)
        seen[row_str] = True

    return data_frame


def load_large_sheet(sheet_obj, in_file, field_list, separator):

    import sys
    reload(sys)
    sys.setdefaultencoding('utf-8')
    import io

    sheet_obj["fields"] = []
    sheet_obj["data"] = []
    seen = {}
    field_ind_list = []
    
    #with io.open(in_file, "r") as FR:
    with io.open(in_file, "r", encoding="utf-8-sig",errors="ignore") as FR:
        rowcount = 0
        ncols = 0
        for line in FR:
            row = line.strip().split(separator)
            #row[0] = row[0].split("\"")[1]
            row[0] = row[0].replace("\"", "")
            row[-1] = row[-1].replace("\"", "")
            rowcount += 1
            if rowcount == 1:
                for j in xrange(0, len(row)):
                    row[j] = row[j].strip().replace("\"", "")
                bad_fields = []
                for f in field_list:
                    if f not in row:
                        bad_fields.append(f)
                if bad_fields != []:
                    print "input file:",  in_file
                    print "fields in file:", row
                    print "fields for extraction:", field_list
                    print "non matching fields: ", bad_fields
                    sys.exit()

                
                #capture number of columns here
                ncols = len(row)
                field_list = row if field_list == [] else field_list
                for f in field_list:
                    field_index = row.index(f)
                    if field_index != -1:
                        sheet_obj["fields"].append(f)
                        field_ind_list.append(field_index)
            else:
                #make sure every row has ncols columns
                if len(row) != ncols:
                    print "bad row %s" % (rowcount)
                    print row
                    sys.exit()
                new_row = []
                for j in field_ind_list:
                    new_row.append(row[j].strip())
                if json.dumps(new_row) not in seen:
                    sheet_obj["data"].append(new_row)
                    seen[json.dumps(new_row)] = True
    return


def load_sheet_as_dict(sheet_obj, in_file, separator, anchor_field):


    seen = {}
    
    if "fields" not in sheet_obj:
        sheet_obj["fields"] = []
    if "data" not in sheet_obj:
        sheet_obj["data"] = {}


    f_list = []
    with open(in_file, 'r') as FR:
        csv_grid = csv.reader(FR, delimiter=separator, quotechar='\"')
        row_count = 0
        for row in csv_grid:
            if json.dumps(row) in seen:
                continue
            seen[json.dumps(row)] = True
            row_count += 1
            if row_count == 1:
                f_list = row
                for j in xrange(0, len(row)):
                    if row[j] == anchor_field:
                        continue
                    sheet_obj["fields"].append(row[j].strip().replace("\"", ""))
            else:
                new_row = []
                for j in xrange(0, len(row)):
                    if f_list[j] == anchor_field:
                        continue
                    new_row.append(row[j].strip().replace("\"", "`"))
                main_id = row[f_list.index(anchor_field)]
                if main_id not in sheet_obj["data"]:
                    sheet_obj["data"][main_id] = [] 
                sheet_obj["data"][main_id].append(new_row)
    return


def load_sheet_from_json(sheet_obj, in_file, field_list):

    sheet_obj["fields"] = field_list
    sheet_obj["data"] = []
    doc = json.loads(open(in_file, "r").read())
    for obj in doc:
        val_dict = {}
        for prop_lineage in field_list:
            val = obj
            for prop in prop_lineage.split("."):
                if prop in val:
                    val = val[prop]
                else:
                    val = ""
            val_dict[prop_lineage] = val
        row = []
        for f in field_list:
            val = str(val_dict[f]) if f in val_dict else ""
            row.append(val)
        sheet_obj["data"].append(row)
    
    return

def load_sheet(sheet_obj, in_file, field_list, separator):

    seen = {}
    sheet_obj["fields"] = []
    sheet_obj["data"] = []
    field_ind_list = []
    with open(in_file, 'r') as FR:
        csv_grid = csv.reader(FR, delimiter=",", quotechar='\"')
        row_count = 0
        ncols = 0
        for row in csv_grid:
            if json.dumps(row) in seen:
                continue
            seen[json.dumps(row)] = True
            row_count += 1
            for j in xrange(0, len(row)):
                row[j] = row[j].replace("\"", "`")
            if row_count == 1:
                ncols = len(row)
                for j in xrange(0, len(row)):
                    f = row[j].strip()
                    if f in field_list:
                        field_ind_list.append(j)
                        sheet_obj["fields"].append(f)
            else:
                #make sure every row has ncols columns
                if len(row) != ncols:
                    continue
                new_row = []
                for j in field_ind_list:
                    new_row.append(row[j].strip())
                sheet_obj["data"].append(new_row)

    
    
    return


def load_workbook(workbook_obj, fileset_objlist, separator):

    for obj in fileset_objlist:
        for file_name in obj["filenamelist"]:
            in_file = obj["dir"] + file_name
            workbook_obj["sheets"][file_name] = {}
            load_sheet(workbook_obj["sheets"][file_name], in_file, ",")

    return



def join_tables(df_one, df_two, anchor_fields):
'''
###
Join Tables
###
###Summary

Combine multiple separate tables based on a given anchor field in each table.

Input:
#######
    * anchorfields : A field from each table that is used to compare then join.

Note:
#######
    * The anchor field must be the first column in each table.

###Format:

{
	"id":#,
	"type":"output",
	"inputtables":[#,#],
	"action":{
		"name":"jointables",
		"anchorfields":["example_field_1","example_field_2"]
		}
	}
}
'''
    df_three = {"fields":[], "data":[]}
    df_three["fields"] = []
    for f in df_one["fields"]:
        df_three["fields"].append(f)

    for f in df_two["fields"]:
        if f != anchor_fields[1]:
            df_three["fields"].append(f)

    row_dict_two = {}
    for row in df_two["data"]:
        anchor_id = row[df_two["fields"].index(anchor_fields[1])]
        if anchor_id not in row_dict_two:
            row_dict_two[anchor_id] = []
        newrow = []
        for j in xrange(0, len(df_two["fields"])):
            if df_two["fields"][j] != anchor_fields[1]:
                newrow.append(row[j])
        row_dict_two[anchor_id].append(newrow)

    row_two_empty = []
    for j in xrange(1, len(df_two["fields"])):
        row_two_empty.append("")


    for row_one in df_one["data"]:
        anchor_id = row_one[df_one["fields"].index(anchor_fields[0])]
        if anchor_id in row_dict_two:
            for row_two in row_dict_two[anchor_id]:
                df_three["data"].append(row_one + row_two)
        else:
            df_three["data"].append(row_one + row_two_empty)


    return df_three







