# Use pdflatex
$pdflatex = 'pdflatex %O %S';

# PDF mode
$pdf_mode = 1;

# Use biber instead of bibtex
$bibtex_use = 2;

# Only keep PDF, delete everything else automatically
@generated_exts = qw(aux bbl bcff blg run.xml out toc lof lot log fls fdb_latexmk nav snm synctex.gz);
$clean_ext = join(' ', @generated_exts);

sub clean_generated_files {
    foreach my $ext (@generated_exts) {
        unlink glob("*.$ext");
    }
}
END { clean_generated_files(); }
