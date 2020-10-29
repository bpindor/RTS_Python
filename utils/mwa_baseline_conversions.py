from numpy import floor
from astropy.io import fits
#global ant2rts, ant2tile,rts2ants,rts2tile
n_ants = 128
n_bl = 0
bl_pairs = []
for i in range(n_ants):
    for j in range(i+1):
        n_bl += 1
        bl_pairs.append((i,j))

rawbl2rawantennas = dict(zip(range(len(bl_pairs)),bl_pairs))    
rawantennas2rawbl = dict(zip(bl_pairs,range(len(bl_pairs))))
# Raw antenna number to rts antenna number

#floor(index/4) + index%4*16 = antenna

raw_ants2rts_ants = [] 
for i in range(int(n_ants/2)):
    raw_ants2rts_ants.append(int(floor(i/4)) + (i%4) * 16)
for i in range(int(n_ants/2),n_ants):
    raw_ants2rts_ants.append(raw_ants2rts_ants[i-int(n_ants/2)] + 64)

rts_ants2raw_ants = dict(zip(raw_ants2rts_ants,range(n_ants)))


# RTS baseline ordering
# rts_baseline = cotterMapping[ant1] + (cotterMapping[ant2] * (cotterMapping[ant2]-1)) / 2;
n_bl = 0
bl_pairs = []
rts_baselines = []
for i in range(n_ants):
    for j in range(i+1,n_ants):
        bl_pairs.append((i,j))
        rts_baselines.append(i + (j * (j-1))/2)

rtsbl2rtsantennas = dict(zip(rts_baselines,bl_pairs))

# Cotter mwaf ordering 00 01 02 ... 0128 11 12 

n_bl = 0
bl_pairs = []
cotter_baselines = []
for i in range(n_ants):
    for j in range(i,n_ants):
        bl_pairs.append((i,j))
cotterbl2ants = dict(zip(range(len(bl_pairs)),bl_pairs))        
ants2cotterbl = dict(zip(bl_pairs,range(len(bl_pairs))))        

def rtsbl2rtsants(rts_bl):
    return rtsbl2rtsantennas[rts_bl]

def rtsants2rtsbl(ant1,ant2):

    if(ant1 < ant2):
        return (ant1 + (ant2 * (ant2-1))/2)
    else:
        return (ant2 + (ant1 * (ant1-1))/2)
        

def raw2rts(raw_bl):

    # calculate the indicies of raw baseline antennas
    # baselines are packed
    # 00 10 11 20 21 22 ...

    raw_ants = rawbl2rawantennas[raw_bl]
    print(raw_ants, raw_ants2rts_ants[raw_ants[0]],raw_ants2rts_ants[raw_ants[1]])
    rts_ants = (raw_ants2rts_ants[raw_ants[0]],raw_ants2rts_ants[raw_ants[1]])
    print(rts_ants)
    rts_baseline = rtsants2rtsbl(rts_ants[0],rts_ants[1])
    print(rts_baseline)
    return rts_baseline
    
    
def rts2raw(rts_bl):

    print(rts_bl)
    rts_ants = rtsbl2rtsantennas[rts_bl]
    print(rts_ants)
    raw_ants = (rts_ants2raw_ants[rts_ants[1]],rts_ants2raw_ants[rts_ants[0]])
    print(raw_ants)
    raw_bl = rawantennas2rawbl[raw_ants]
    print(raw_bl)
    return raw_bl
                     
def loadCotterMapping(metafile):

    meta_file = fits.open(metafile)
    antenna_n = meta_file[1].data['Antenna']
    tile_n = meta_file[1].data['Tile']

    global ant2rts 
    ants2rts = dict(zip(antenna_n[::2],range(n_ants)))
    global rts2ants
    rts2ants = dict(zip(range(n_ants),antenna_n[::2]))
    global ant2tile 
    ant2tile = dict(zip(antenna_n[::2],tile_n[::2]))
    global tile2ant
    tile2ant = dict(zip(tile_n[::2],antenna_n[::2]))
    global rts2tile
    rts2tile = dict(zip(range(n_ants),tile_n[::2]))
    global tile2rts_d
    tile2rts_d = dict(zip(tile_n[::2],range(n_ants)))

def cotterbl2tiles(cotter_bl):
    cotter_ants = cotterbl2ants[cotter_bl]
    tiles = (ant2tile[cotter_ants[0]],ant2tile[cotter_ants[1]])
    return tiles

def tiles2Cotterbl(tiles):
    cotter_ants = (tile2ant[tiles[0]],tile2ant[tiles[1]])
    cotter_bl = ants2cotterbl[cotter_ants]
    return cotter_bl

def rtsbl2Cotter(rts_bl):

    rts_ants = rtsbl2rtsantennas[rts_bl]
    cotter_ants = (rts2ants[rts_ants[0]],rts2ants[rts_ants[1]])
    if(cotter_ants[0] < cotter_ants[1]):
        cotter_bl = ants2cotterbl[cotter_ants]
    else:
        cotter_bl = ants2cotterbl[cotter_ants[::-1]]
    return cotter_bl
    
def raw2tiles(raw_bl):

    rts_bl = raw2rts(raw_bl)
    rts_ants = rtsbl2rtsantennas[rts_bl]
    tiles = (rts2tile[rts_ants[0]],rts2tile[rts_ants[1]])
    print(tiles)
    return tiles

def tile2rts(ant):

    return tile2rts_d[ant] 

def rtsbaselinesForAnt(ant1):
    
    all_baselines = []

    for ant2 in range(n_ants):
        if(ant1 != ant2):
            all_baselines.append(rtsants2rtsbl(ant1,ant2))

    return all_baselines


    
    
    


    




    
    
