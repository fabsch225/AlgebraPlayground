\\ ecpp_cert.gp -- produce and export an ECPP primality certificate in PARI/GP
\\ Save this file and run: gp -q ecpp_cert.gp

n = 2^255 - 19;
print("n = ", n);

cert = primecert(n);

ok = primecertisvalid(cert);        \\ returns 1 for valid, 0 otherwise
print("primecertisvalid = ", ok);

s_human = primecertexport(cert, 0);    \\ human-readable ECPP export
s_primo = primecertexport(cert, 1);    \\ Primo format (if needed)

write("2p25519_ecpp_cert.txt", s_human);
write("2p25519_ecpp_cert_primo.txt", s_primo);

print("Wrote: 2p25519_ecpp_cert.txt and 2p25519_ecpp_cert_primo.txt");
