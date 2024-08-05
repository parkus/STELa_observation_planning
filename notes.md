# 2024-08-04

## gaia vs G mag
supposedly these are the same since SIMBAD refs gaia EDR3 as the source of its G mags, but I noticed there are differences between the sy_gaiamag and sy_gmag columns in the exoplanets table
I'll use the latest values from SIMBAD

## time constraint for downselection
Question I drafted then decided to ask later (after I see if a time constraint is problematic)
Separate-ish question 3: Our downselection *between* (rather than within) cycles also needs to be considered. In the proposal, we included as a special requirement that our observations be completed with sufficient lead time ahead of the next cycle's phase II deadline to enable the downselection step. However, alternatives might be possible, like delaying our phase II submission, dividing our allocation between the regular and mid-cycle deadlines, or allowing some of our allocations to bleed over into subsequent cycles (as I suspect will happen anyway due to failed orbits). Would one of the alternates be preferable to STScI or should we proceed with a time constraint to enable downselection lead time as stated in the proposal? 

## exposures
need put spectra for brightest of each spectral type into ETC to verify no overbrightness concerns
have to do ETC for each spectral type to get ACQ exposure times... or do I? why not interpolate over a table of Teff vs G. could do this to test for overbrightness concerns too. We could also do the worst cases for early and mid M stars

For ACQ I need exposure time and to ensure no overbrightness. I guess that will come with exposure times. If they get too low, I will switch to a different filter.

For science, I need *buffer time* predictions along with checking for no overbrightness. Perhaps I can make a grid of UV activity level, actually maybe GALEX UV which I can then either retrieve or estimate. From that I can estimate buffer times and determine where overbrightness concerns happen.  

I think I should first see if there are any targets where buffer times or overbrightness is a concern. If so, I can write the machinery to create more detailed spectra for those specific targets. 

Running the brightest FUV spectra through G140L ETC now...
F2 star at FUV=14.7 no violation
K1 at FUV=19.3 no violation
so don't have to worry about bright violations here

Now trying G140M
F2 star at FUV=14.7 no violation
K1 at FUV=19.3 no violation
M4 at FUV=20 no violation

Yay! so I just need to compute buffer times



# 2024-08-03

## bad initial target selection!?!?! no it was okay
looks like the mask I created to do the first target cut was itself a masked array, which resulted in unpredictable behavior
Masked columns resulted finding the TOIs Shreryas identified as false positives and the flag_snr_okay column
redoing the same cuts and being careful of how masking is handled, I get a list of 122 targets. same as before. phew.

It looks like there were some great targets that were cut from the initial set and I do not understand why. Here are the top 10 from a broader cut ranked by the optimistic SNR estimate from ethan. 

hostname st_spectype pl_rade  Cnom Copt Enom  Eopt cut_from_proposal
                     earthRad                                       
-------- ----------- -------- ---- ---- ---- ----- -----------------
HD 95338      K0.5 V    3.890 20.1 39.3 90.0 123.0              True
TOI-2194          --    1.989 17.8 34.7 83.2 117.8              True
  GJ 143        K4.5    2.610  5.5 38.6 70.6 106.8             False
      --          --    2.537 36.6 69.5 41.5  74.0             False
      --          --    2.964  6.1 20.3 28.5  64.8             False
  GJ 357      M2.5 V    1.200 11.4 17.3 25.8  34.6             False
HD 63935          --    2.900 10.7 28.7 11.7  30.6              True
      --          --    2.629  4.1 13.7 14.1  26.2             False
  GJ 143        K4.5    0.892  1.3  8.5 16.5  24.9              True
  K2-415        M5 V    1.015  9.3 23.7 15.9  22.6              True
(E = ethan, C = crude, opt = optimistic, nom = nominal)

For most, it was the ratio of Ethan's SNR estimate to the crude estimate. It looks like it has to do with the existence of some duplicate(?) columns: 'lya_transit_snr_sim' vs 'pl_lya_tnst_snr_optimistic'
Where do these come from?
'pl_lya_tnst_snr' -- crude estimate. I left a note by it saying that it way overestimates the SNR. This is the problem.
'lya_transit_snr_sim' -- these are the actual simulations from ethan 

Redoing this, things make a lot more sense. There aren't any targets that have an optimistic simulated SNR that is really impressive. The best is 3.6. There are quite a few targets where my crude SNR estimate far exceeds ethan's, but this could likely indicate that the outflow is highly irradiated so that I am vastly overstimating the occulter size. Still, I think these make interesting targets. I'm going to replace the A, WD, and F stars with the 6 planets where the crude SNR estimate is >10x the simulation and planet is above the gap + TOI-431 since the optimistic simulated SNR is > 3. 

## bad types?
There is an A star in the sample too. Better get rid of anything more massive than F or at least consider doing so. There is only the one though. Actually there are a few objects where simbad and my table have much different spectral types, or both give a type we should consider removing. Looking into them:

√ TOI-2056 A star? has only TIC data, my spectype of F3 based on that, Simbad's A type is graded worst quality, so let's stick with TIC
√ TOI-880 evolved? simbad gives K2 III but graded worst quality, TIC gives MS logg and planets are included in a pub about best candidates for JWST obs, so let's keep
X WASP-189 A star? simbad spt worst quality, but exo archive gives Teff of 8000 K which corresponds to A6 in mamajek table. Using HD 14943 as a proxy (about 50% closer), there will be no flux at Lya. Cull. 
√ EPIC205530323 simbad gives M4 quality D (second to worst), TIC gives a Teff of 3142 but somehow I get a type of A6. dunno what that is about. redoing with the TIC Teff I get M4.5. We should add this one to the M dwarf list. And keep it. 
X WD1856+534 obviously a white dwarf, probably has massive Lya absorption 
√ HD235088 archive gives K2, simbad gives G5 quality E. K2 is the right type. 
~ TOI-2112 A star? A5 quality E in simbad. mamajek table gives F0 based on TIC temp of ~7200 K. Hmmm. I don't love this one. 

Looking at some other F0-F2 stars with HST data
HD432 strong Lya but it is a giant
HD80404 absorption but supergiant
HD156389 has emission, but F2 type is unreliable
and that's all that is in the archive. hmm. 
HD28568 has data in starcat, listed as F2V, no spectral type in simbad, has clear Lya emission

Checked for stars where simbad gives an M type but I didn't have M in my table. EPIC star is the only one. 

## early Fs
TOI-2056 F3 low SNR predicted. 1 planet neptune 10 d. 
TOI-5099 F0 1 planet, jupiter, 14 days. 90 pc. 
TOI-5807 F2 
TOI-2112 F0 2 planets, neptunes, at 14 and 155 days. 86 pc. 

keep for now but make these the first to cut

# 2024-08-01
I noticed for HD 15337 (target of Tail Length proposal), it is in the phase I list even though it has observations, probably because it does not have a G140M observation. But it does have G140L, so we could reallocate that orbit. 

In general, I should probably search again to find stars that have G140M or G140L data in the archive/planned. Ah, but I must be careful because the planned data might be proprietary for some period of time, so some by-hand checking might be key. Then I should also search the MAST abstracts for "lyman". Doing this just now turns up 19 programs, but the only stellar ones are the ones our group proposed, with the exception of one for old G stars by Brian Wood. Searching "lya" turns up 6 programs, none about stars. Ly-a turns up 3 more, still no stars. 

Swapping targets is a minor change, so that's great! 


# 2024-07
see notes in the target search possible additions folder for that info