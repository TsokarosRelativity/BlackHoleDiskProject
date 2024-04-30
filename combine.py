import pandas as pd
import time 

print(f'STARTING')
s1 = time.time()
df_nh = pd.read_hdf('3D_data/not_horizon_data.h5', key='df')
df_h = pd.read_hdf('3D_data/horizon_data.h5', key='df')
s2 = time.time()
print(f'loaded in data in {s2 - s1}sec')

ret_df = pd.concat([df_nh, df_h], axis=0, ignore_index=True)
s3 = time.time()
print(f'concatenated in {s3-s2}sec')

ret_df.to_hdf('3D_data/all_data.h5', key='df', mode='w')
print(f'saved in {time.time()-s3}sec')
