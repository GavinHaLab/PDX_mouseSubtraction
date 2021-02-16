#!/usr/bin/env python
# coding: utf-8

# In[1]:


#define references

human_ref_path = '/fh/fast/ha_g/grp/reference/GRCh38/GRCh38.fa'

mouse_ref_path = '/fh/fast/ha_g/grp/reference/Mus_musculus_NCBI_GRCm38/NCBI/GRCm38/Sequence/WholeGenomeFasta/genome.fa'

out_file = '../fasta/GRCh38_plus_NCBI_GRCm38.fa'


# In[2]:


#read in the files
with open(human_ref_path,'r') as f:
    human_ref = f.readlines()
    
with open(mouse_ref_path,'r') as f:
    mouse_ref = f.readlines()


# In[3]:


#view the first few lines
print('# of lines in human:',len(human_ref))
print(human_ref[0:10])
print('\n')
print('# of lines in mouse:',len(mouse_ref))
print(mouse_ref[0:10])


# In[4]:


#rename the mouse contigs
for i in range(len(mouse_ref)):
    if mouse_ref[i].startswith('>'):
        mouse_ref[i] = mouse_ref[i].strip('\n')+('_NCBI_GRCm38\n')


# In[5]:


#export to new file
with open (out_file, 'w+') as f:
    f.write(''.join(human_ref))
    f.write(''.join(mouse_ref))


# In[ ]:




