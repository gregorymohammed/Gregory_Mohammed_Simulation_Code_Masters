# =====================================================
# Gregory Mohammed
# MSc. Physics
# README file for NetBeans Projects folder.
# =====================================================
# 1st April, 2014
#
# Godunov-v.6.1.6 is a copy of Godunov-v.6.1.5. It 
# works in exactly the same way except it is using
# mks units for the KKT data.
# =====================================================
# 27th March, 2014
#
# Godunov-v.6.1.3 does not work.
# Godunov-v.6.1.4 works perfectly.
# Godunov-v.6.1.5 works perfectly. Use this to copy from
# for doing further tests.
# =====================================================
# 2nd September, 2013
#
# Copy of Godunov-v.6.0.3 is Godunov-v.6.1.0. This
# version deals with cleaning up remaining bugs, which
# are to do with usage. The actual computation is 
# correct.
# =====================================================
# 12th August, 2013
#
# My birthday! Godunov-v.6.0.3 is the version which now
# works perfectly. Godunov-v.6.0.4 is demoted to old
# version, which functions incorrectly.
# =====================================================
# 8th August, 2013
#
# Fixed bugs in Godunov-v.6.0.3 using the latest version,
# Godunov-v.6.0.4!! No blips anymore, and standard 
# deviations are correct! Now for some serious 
# experimentation!
# =====================================================
# 21st November, 2012
#
# Godunov-v.6.0.2 contains methods for reading all the .csv
# files in the working directory and picking out only
# the files with approximate solutions. These names
# are then written to a file which is used for gnuplot
# to do an animation. GIMP is used to merge the animation
# into a .gif
# =====================================================
# 30th September, 2012
# 
# BUG! In 6.0.1 - You have to run a simulation with 
# showExact = 1 first then again with showExact = 0,
# ONLY if you do not want to see the exact solution.
# This is only relevant for Sod data, where the 
# analytic solution exists.
# =====================================================
# 31 August, 2012
# 
# Godunov-v.6.0.1 => This works perfectly. Very detailed 
# param.dat file with numerous comments. Numerous switches.
# Neutrinos properly implemented. 
#
# This involved corrections to the neutrino model using
# Janka's work. The correction was in the kappa used.
# Had to implement this in this code, so do not use
# previous versions.
# =====================================================
# 21 August, 2012
#
# Godunov-v.6.0.0 => Extensive work done. All the 
# problems with the previous version are solved.
# Had to make numerous modifications and tests to
# make the code comply with the Kuroda et al. data.
# Program works perfectly. The param.dat file is 
# extensive, and needs to be understood before using.
# =====================================================
# 20 August, 2012
#
# Godunov-v.5.0.6 => Eureka! Code works perfectly. Great
# simulation with Kuroda et al.'s data. Sod Tube also 
# still works. Geometrizations are all correct. Copied
# to Godunov-v.6.0.0. However, neutrino switch isn't 
# working.
# =====================================================
# 15 August, 2012
#
# Godunov-v.5.0.6 => 
#
# This version incorporates the Kuroda et al. data,
# which is used to make a KKT set of initial conditions.
# It also uses a geometrizer class to geometrize the
# input parameters, and also to un-geometrize the
# output so that the graphs show data in SI units.
#
# Is 5.0.5, with a gnuplot scripter
# which generates gnuplot batch files for use with
# gnuplot.
# =====================================================
# 4 August, 2012
#
# Godunov-v.5.0.5 => Is version 5.0.4, with a tunable
# source term with a range of 0 < sourceTuner <= 1.
# =====================================================
# 4 August, 2012
#
# Godunov-v-5.0.4 => Exact Riemann solver is implemented,
# that is, to generate the analytic output to overlay on
# the approximate output. Standard deviation is also
# implemented. Have to put tweaking factor for source
# terms. Have to set up gnuplot plotting.
# =====================================================
# 2 August, 2012
#
# Godunov-v-5.0.3 => Clean code, corrected errors in 
# usage of class eqs and IO. Runtime better by 700%.
# =====================================================
# 30 July, 2012
# Godunov-v-5.0.2 => This is fully functional. Use this
# code for testing and results.
# =====================================================
# 28 July, 2012
# 
# Godunov-v-5.0.2 => Copy of 5.0.1 for performing
# major work with baseVariables, see below.
#
#  - 29 July, 2012
# Discovered several errors. Fixed. Implemented PJM's 
# high resolution methods.
# =====================================================
# 27 July, 2012
#
# Godunov-v-5.0.1 => Modified 5.0.0 by PJM.
# 
# Major problem with this one is that the cur and old
# do not communicate across classes. I think this means 
# that cur and old have to be passed around, rather than
# new'ed in a cell. This version is kept as is for 
# replication. 
# =====================================================
# 12 July, 2012
#
# Godunov-v-5.0.0 => This version is Godunov-v-4.0.5
# with the neutrino class.
# Points to Note:
# 1. Data files use spaces as separators for use with
#    gnuplot.
# 2. Source terms are mulitplied by delta t.
# =====================================================
# 12 July, 2012
#
# Godunov-v-4.0.5 => Is the latest and works perfectly.
# This version has a timer called StopWatch, and writes
# the runtime in milliseconds to the console. This 
# version does not have the neutrino equations.
#
# =====================================================
# 17 March, 2012
#
# Godunov-v-2.0.4 => Works for Upwind, implements MM,
# but does not work for MM.
#
# Godunov-v-2.0.5 => Works but MM does not solve.
#
# Godunov-v-3.0.0 => Junk.
#
# =====================================================
# 4 March, 2012
#
# The project Godunov-v-2.0.4 works and implements
# the Marti-Muller Riemann solver. But, the MM does
# not work!!
# I have created a new project called:
# -----------------------------------------------------
# FiniteVolumeGodunovUpwindCode-v-1.0.0
# -----------------------------------------------------
# which is Godunov-v-2.0.5. This seperates and preserves
# the working code which implements Marti-Muller from
# the one which implements Rezzolla-Zanotti.
# =====================================================
# 12 April, 2011
#
# =====================================================
# The project Godunov-v.2.0.3 works for the Upwinding
# method on a smoothed shock and a real shock.
# -----------------------------------------------------
# The Godunov part, where the Riemann solver is
# implemented, is not proven to work in this 
# version.
# -----------------------------------------------------
# The project Godunov-v.2.0.4 works for both the
# Upwinding and the Jump, for both the smoothed and
# un-smoothed initial conditions.
# =====================================================
# NOTE: Godunov-v.3.0.0 does not work! 
# Ignore this version (3.0.0); a futile exercise.
# =====================================================
# Next major code changes will start at 4.0.0.
# =====================================================
