Summary
-------------------------------------------------------------------------------------
Recipe files are used to describe the actions that the script "make-datasets.py" will perform on a ".csv" dataset file, to process the dataset for integration into a OncoMX or Glygen. 

Recipe files are in ".json" format. For more information on the ".json" format, go to www.json.org

The recipe consists of two major sections: "outputfiles" and "tablelist". 

Put simply, in "tablelist" you will create tables from the list of "actions" that can be performed. 

You will then use "outputfiles" to denote what tables you wish to export and the content to include from those tables. 

The "outputfile" is typically the reformatted dataset ".csv" that will then move forward in the pipeline to eventually be included on the frontend of OncoMX or Glygen. 
