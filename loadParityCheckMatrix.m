function H=loadParityCheckMatrix(enc_input_word)
% Function: loadParityCheckMatrix
%
% H=loadParityCheckMatrix(enc_input_word)
%
% Author: Jos� Carlos Aradillas 
%
% Contact: murillo@us.es
%
% Created 27/06/2016
%
% Description: This function loads a parity check matrix for the LDPC
% encoder and decoder depending on input parameter enc_input_word

switch enc_input_word
    case 128,
        load ../ParityCheckMatrix/H128r256;
    case 256,
        load ../ParityCheckMatrix/H256r512;
    case 512,
        load ../ParityCheckMatrix/H512r1024;
    case 1024,
        load ../ParityCheckMatrix/H1024r2048;
    case 2048,
        load ../ParityCheckMatrix/H2048r4096;
    case 4096,
        load ../ParityCheckMatrix/H4096r8192;
    otherwise
        display('Invalid ParityCheckMatrix')
        return
end
end