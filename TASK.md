# BACKGROUND
- I am designing a plasmid vector for my biology course 
- I want to predict which promoter sequence will be best for a plasmid vector inserted into chassis C. acnes
- I have a promoter, P3 from S. aures, but it is a weak promoter, and so I want to alter the -10 and -35 regions (this has been done in a paper previously: https://doi.org/10.1128/spectrum.01829-23)
- Altering the promoter regions in these 2 spots increased the GFP activity by the plasmid 
- I came up with a list of strong promoters found in C. acnes, and to predict which will be the strongest I was told by my professor to use LLMs. 
 
# PROBLEMS TO SOLVE
- Was told to use hugging face to find some pre-trained LLMs or train them myself. I found one for promoters (https://huggingface.co/blog/hugging-science/promoter-gpt#if-dna-is-truly-a-language-then-we-should-be-able-to-teach-transformers-how-to-write-it). I want to train it off of either the genome of C. acnes (maybe? idk? I have no idea what is going on) or the promoter regions? 
- I think I just want to give the LLM my promoter regions, and then ask what promoter regions will be the best in C. acnes
- The model I found has a specific data set and I don't know how I am supposed to alter the data to train my model since I don't really understand what all the columns are. So basically I found some code but I don't know how to edit the code or edit the data so that I can feed it my data. 

Based on the background and problems to solve can you recommend some next steps.
