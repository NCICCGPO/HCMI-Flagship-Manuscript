import os
import sys
import argparse
import numpy as np 
import pandas as pd 


#handle arguments
def parse_args():
    parser = argparse.ArgumentParser(description="Compute consensus purity/ploidy")
    parser.add_argument("--input", help="Input file")
    parser.add_argument("--output_path" , help="Output path")
    args = parser.parse_args()
    return args


#This function is useful for eliminating callers that defer
def is_numeric(value):
    try:
        float(value)
        return True
    except ValueError:
        return False


#remove all non-numeric values
def convert_list_to_floats(lst):
    return [float(x) for x in lst if str(x).replace('.', '', 1).isdigit()]


#remove the \n from the output dataframe
def remove_newline(df):
    df = df.replace('\n','', regex=True)
    return df
    

#Define a function to test if all purity callers fall within the consensus window
#This function will return a consensus purity value if 
    #A) Two or fewer callrs fall outside the window
    #B) At least four callers fall within the window
def test_all_fall_within_purity_window(purity_values):

    #Stratify purity_values into two lists, one with values in range, the other with values out of range. If there are fewer than 4 values in range, delete the furthest value in purity_values and re-test. If purity_values reaches 4 and there are still values out of range, then return NaN
    while len(purity_values) >= 3:
        mean_purity = sum(purity_values)/len(purity_values) 
        consensus_range = (mean_purity - 0.2, mean_purity + 0.2) #Get the range for purity consensus. All callers must fall within this window.
        in_range = [] 
        out_range = [] 
        lower_bound, upper_bound = consensus_range 
        for val in purity_values: 
            if lower_bound <= val <= upper_bound:
                in_range.append(val) 
            else: 
                out_range.append(val) 
        if len(in_range) >= 3: #If we have enough callers to be confident in the call, then return a consensus value
            consensus_purity_value = sum(in_range)/len(in_range) 
            return consensus_purity_value 
        elif len(in_range) < 3: #If we don't have enough callers to be confident in the call, then delete the value in purity_values that is furthest from mean_purity, then re-build in_range and out_range
            furthest_value = min(purity_values, key=lambda x:abs(x-mean_purity))
            purity_values.remove(furthest_value)
            
    if len(purity_values) < 3:
        return np.nan




#Define a function to classify ploidy values and compute consensus ploidy
def compute_consensus_ploidy(ploidy_values, ploidy_absolute, ploidy_purple):
    
    #First, classify the values into haploid, diploid, triploid, and tetraploid (which also includes > tetraploid)
    classified_ploidy = []
    for ploidy_value in ploidy_values:
        if ploidy_value < 1.5:
            classified_ploidy.append("haploid")
        elif 1.5 <= ploidy_value < 2.5:
            classified_ploidy.append("diploid")
        elif 2.5 <= ploidy_value < 3.5:
            classified_ploidy.append("triploid")
        elif ploidy_value >= 3.5:
            classified_ploidy.append("tetraploid")
        else:
            classified_ploidy.append(np.nan)
            

    #Count how many times each classification appears in the list, if more than three times then pick that classification
    if classified_ploidy.count('haploid') >= 3:
        filtered_ploidy_values = [value for value in ploidy_values if value < 1.5] 
        mean_ploidy = np.mean(filtered_ploidy_values)
    elif classified_ploidy.count('diploid') >= 3:
        filtered_ploidy_values = [value for value in ploidy_values if 1.5 <= value < 2.5]
        mean_ploidy = np.mean(filtered_ploidy_values)
    elif classified_ploidy.count('triploid') >= 3:
        filtered_ploidy_values = [value for value in ploidy_values if 2.5 <= value < 3.5]
        mean_ploidy = np.mean(filtered_ploidy_values)
    elif classified_ploidy.count('tetraploid') >= 3:
        filtered_ploidy_values = [value for value in ploidy_values if value >= 3.5] 
        mean_ploidy = np.mean(filtered_ploidy_values)
    else:
        mean_ploidy = np.nan
    
    
    ###If no classification appears more than twice, then pick absolute. If no value for absolute, then pick purple. If no value for purple, then set ploidy as NaN
    ploidy_class_occurrences = {}
    for ploidy_class in classified_ploidy:
        if ploidy_class in ploidy_class_occurrences:
            ploidy_class_occurrences[ploidy_class] += 1
        else:
            ploidy_class_occurrences[ploidy_class] = 1
            
    any_class_occurs_more_than_twice = any(count > 2 for count in ploidy_class_occurrences.values())
    
    if any_class_occurs_more_than_twice == False:
        if is_numeric(ploidy_absolute):
            mean_ploidy = ploidy_absolute
        elif not is_numeric(ploidy_absolute) and is_numeric(ploidy_purple):
            mean_ploidy = ploidy_purple
        else:
            mean_ploidy = np.nan
            
    return mean_ploidy




#loop through all lines in the input file
def load_and_process_df(input):
    with open(input, 'r') as file:

        #Get a list to store all consensus purity values
        final_consensus_purity_list = []
        final_consensus_ploidy_list = []
        
        #Create lists to store all of the information for the output pandas df
        id3_list = []
        sample_id_list = []
        sample_type_list = []
        purity_ascat_list = []
        purity_absolute_list = []
        purity_remixt_list = []
        purity_purple_list = []
        purity_estimate_list = []
        purity_hatchet_list = []
        ploidy_ascat_list = []
        ploidy_absolute_list = []
        ploidy_remixt_list = []
        ploidy_purple_list = []
        ploidy_hatchet_list = []
        

        for line in file:
            
            split_line=line.split('\t')

            if split_line[0] == "id3": #Skip the header
                continue
            else:
                id3=split_line[0]
                sample_id=split_line[1]
                sample_type=split_line[2]
                purity_ascat=split_line[3]
                purity_absolute=split_line[4]
                purity_remixt=split_line[5]
                purity_purple=split_line[6]
                purity_estimate=split_line[7]
                purity_hatchet=split_line[8].replace('\n', '') #Remove the newline character
                ploidy_ascat=split_line[9]
                ploidy_absolute=split_line[10]
                ploidy_remixt=split_line[11]
                ploidy_purple=split_line[12]
                ploidy_hatchet=split_line[13].replace('\n', '') #Remove the newline character

                #Get a list with all purity and ploidy values
                purity_list=[purity_ascat, purity_absolute, purity_remixt, purity_purple, purity_estimate, purity_hatchet]
                ploidy_list=[ploidy_ascat, ploidy_absolute, ploidy_remixt, ploidy_purple, ploidy_hatchet]
                purity_list_numeric = convert_list_to_floats(purity_list)
                ploidy_list_numeric = convert_list_to_floats(ploidy_list)
                
                #First, test if all callers fall within the window, if so, then purity is set.
                consensus_purity = test_all_fall_within_purity_window(purity_list_numeric)

                #If purity is not set, then pick absolute as the purity value. If there is no absolute value then pick purple
                if np.isnan(consensus_purity):
                    if is_numeric(purity_absolute):
                        consensus_purity = purity_absolute
                    elif not is_numeric(purity_absolute) and is_numeric(purity_purple):
                        consensus_purity = purity_purple
                    else:
                        consensus_purity = np.nan
                
                #Test if all callers fall within the ploidy window
                consensus_ploidy = compute_consensus_ploidy(ploidy_values = ploidy_list_numeric, ploidy_absolute = ploidy_absolute, ploidy_purple = ploidy_purple)
                
                #append to the lists outside of the loop
                final_consensus_purity_list.append(consensus_purity)
                final_consensus_ploidy_list.append(consensus_ploidy)
                
                #Append the other values to the lists
                id3_list.append(id3)
                sample_id_list.append(sample_id)
                sample_type_list.append(sample_type)
                purity_ascat_list.append(purity_ascat)
                purity_absolute_list.append(purity_absolute)
                purity_remixt_list.append(purity_remixt)
                purity_purple_list.append(purity_purple)
                purity_estimate_list.append(purity_estimate)
                purity_hatchet_list.append(purity_hatchet)
                ploidy_ascat_list.append(ploidy_ascat)
                ploidy_absolute_list.append(ploidy_absolute)
                ploidy_remixt_list.append(ploidy_remixt)
                ploidy_purple_list.append(ploidy_purple)
                ploidy_hatchet_list.append(ploidy_hatchet)
                
                
        
        #Bind everything together into a single dataframe
        output_df = pd.DataFrame({
            "id3": id3_list,
            "sample_id": sample_id_list,
            "sample_type": sample_type_list,
            "purity_ascat": purity_ascat_list,
            "purity_absolute": purity_absolute_list,
            "purity_remixt": purity_remixt_list,
            "purity_purple": purity_purple_list,
            "purity_estimate": purity_estimate_list,
            "purity_hatchet": purity_hatchet_list,
            "ploidy_ascat": ploidy_ascat_list,
            "ploidy_absolute": ploidy_absolute_list,
            "ploidy_remixt": ploidy_remixt_list,
            "ploidy_purple": ploidy_purple_list,
            "ploidy_hatchet": ploidy_hatchet_list,
            "consensus_purity": final_consensus_purity_list,
            "consensus_ploidy": final_consensus_ploidy_list
        })
               
        #Return everything to outside the loop
        return output_df


#Define a manual curation function
def manual_curation(output_df):
    
    #Make a list of all samples we need to manually curate
    sample_list = ['HCM-CSHL-0907-C50-06A',
                   'HCM-BROD-0019-C25-01B',
                   'HCM-BROD-0578-C15-06A',
                   'HCM-CSHL-0191-C25-01A',
                   'HCM-BROD-0869-C43-85A',
                   'HCM-BROD-0053-C49-85M',
                   'HCM-CSHL-0246-C19-01B',
                   'HCM-BROD-0231-C25-85N',
                   'HCM-CSHL-0082-C25-85A',
                   'HCM-WCMC-0949-C67-85A'
                   'HCM-BROD-0051-C64-06B',
                   'HCM-BROD-0209-C71-85A',
                   'HCM-CSHL-0385-C18-85A',
                   'HCM-CSHL-0075-C25-85A']
    
    #Loop through them and manually set them to the absolute value
    for sample in sample_list:
        output_df.loc[output_df['sample_id'] == sample, 'consensus_purity'] = output_df['purity_absolute']
        output_df.loc[output_df['sample_id'] == sample, 'consensus_ploidy'] = output_df['ploidy_absolute']
    
    return output_df
    


#Define a function to write the output to a file
def write_output(output_df, output_path):
    output_df.to_csv(output_path + "/final_consensus_pp.txt", sep='\t', index=False)     
           


#Define the main function
def main():
    #Load the arguments
    args = parse_args()
    
    #Load the input file and process it
    consensus_pp = load_and_process_df(args.input)
    consensus_pp = remove_newline(consensus_pp)
    consensus_pp = manual_curation(consensus_pp)
    
    #Write the output to a file
    write_output(consensus_pp, args.output_path)
    
  

#Do the thing
if __name__ == "__main__":
    main()










