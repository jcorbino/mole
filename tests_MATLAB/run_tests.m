% Main script to run all the tests

clc
clear;
close all;

disp('Running: Nullity test of Divergence operator...');
run('test1.m');

disp('Running: Nullity test of Gradient operator...');
run('test2.m');

disp('Running: Nullity test of Laplacian operator...');
run('test3.m');

disp('Running: Energy test (Schrödinger equation)...');
run('test4.m');

disp('Running: Elliptic test (Poisson equation)...');
run('test5.m');