Recipe Actions
==============

.. toctree::
    :maxdepth: 3
    :caption: Full list of actions:
    union
    addconstantfield
    makecombofield
    filterin
    filterout
    expandrows
    transposecols
    addnormalizedcol
    addaveragecol
    splitcol
    jointables

union 
------
Summary
^^^^^^^^
Merge multiple tables and remove duplicate lines.

Format
^^^^^^^
Output Format::

    {
        "id":#,
        "type":"output",
        "inputtables":[#,#,#],
        "action":{
            "name":"union"
        }
    }

Notes:
^^^^^^
This action can take multiple inputtable tables. 

See use-case recipe example_part1.json for usage example.


addconstantfield 
-----------------

Summary
^^^^^^^

Add a new field to the table that has the same constant entry for each line.

Format
^^^^^^
Output Format::

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

Notes: 
^^^^^^
This action can take multiple inputtable tables. 

See use-case recipe example_part2.json for usage example.

makecombofield
--------------
Summary
^^^^^^^

Combine two fields into one field. 

Format
^^^^^^
Output Format::

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

Notes:
^^^^^^
This action can take a single inputtable table. 

The above example makes the combined field example_combo_field from example_field_1 and example_field_2. 

The entries in each of the fields will merge into the combined field and be separated by the merge character, in the example the charcater is |.


See use-case recipe example_part1.json for usage example.

filterin
---------

Summary
^^^^^^^

Filter only certain lines from one table into another.

Needs review from Robel.

Format
^^^^^^
Output Format::

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

Notes: 
^^^^^^
This action can take a single inputtable table. 

In the above example, only the lines with the values ABC, -, or an empty field will be included in to a new table.

See use-case recipe example_part2.json for usage example.

Question for Robel:
What is the purpose of "operation"?

filterout
---------
Summary
^^^^^^^

Filter out only certain lines from a table to create a table without those entries.

Needs review from Robel.

Format
^^^^^^
Output Format::

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

Notes: 
^^^^^^
This action can take a single inputtable table. 

In the above example, lines with the values DEF, 123, or 456 in the field example_field_1 will be excluded from a new table, and all other lines with other values for that field will be included.

See use-case recipe example_part1.json for usage example.

Question for Robel:
What is the purpose of "operation"?


expandrows
----------

Summary
^^^^^^^

For a single line that has a field with multiple data points in one field, create multiple lines where each line has a single data point from that field. 

Format
^^^^^^
Output Format::

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

Notes: 
^^^^^^
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

transposecols
-------------

Summary
^^^^^^^

The rows switch with columns and vice versa(?)

Needs review from Robel.

Format
^^^^^^
Output Format::

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

Notes: 
^^^^^^
This action can take a single inputtable table.

See use-case recipe example_part2.json for usage example.

Question for Robel:
How is this action used?

addnormalizedcol
----------------
Summary
^^^^^^^

Create a new column that takes the numerical data from an existing column and normalizes the numbers according to (?)

Needs review from Robel.

Format
^^^^^^    
Output Format::

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

Notes: 
^^^^^^
This action can take a single inputtable table.

See use-case recipe example_part2.json for usage example.

Question for Robel:
How is normalization calculated?


addaveragecol
-------------

Summary
^^^^^^^

Create a new column that takes the numerical data from an existing column and averages the numbers according to (?)

Needs review from Robel.

Format
^^^^^^  
Output Format::
        
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

Notes:
^^^^^^
This action can take a single inputtable table.

See use-case recipe example_part2.json for usage example.

Question for Robel:
How is the average being calculated?



splitcol
--------

Summary
^^^^^^^

Split a single column into two separate columns after specifying a deliminating character that separates the values. 

Format
^^^^^^
Output Format::

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

Notes: 
^^^^^^
This action can take a single inputtable table.

See use-case recipe example_part1.json for usage example.



jointables
----------
Summary
^^^^^^^

Create a new table that uses an anchor field to combine columns from multiple tables. Data with the same anchor field value will be added to the same row in the new field. 

Format
^^^^^^
Output Format::

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

Notes:
^^^^^^
This action can take multiple inputtable table.

See use-case recipe example_part1.json for usage example.
