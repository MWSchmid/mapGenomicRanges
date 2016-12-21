#include <QtCore/QCoreApplication>
#include <QtCore>
#include "../Rcount/source_V2/p502_SOURCE/dataStructure/databaseitem.h"
#include "../Rcount/source_V2/p502_SOURCE/dataStructure/database.h"
#include "../Rcount/source_V2/p502_SOURCE/dataStructure/mappingtreeitem.h"
#include <iostream>
#include <getopt.h>
#include <QTime>

//! genomic region class
// note: we create just once, use setInit to write the info from the file, send it through annotation mapping, then append it to the vector
class region {
public:
    // set at beginning
    QString chrom;
    uint start;
    uint end;
    QString strand;
    QStringList otherData;
    // set in getNucMapping
    QList<databaseItem*> bestmappingfeatures; //holds the pointers to the data vectors

    region(): chrom(""), start(0), end(0), strand(".") { }

    void setInit(QString &chrom, QString &strand, QString &start, QString &end) {
        this->chrom = chrom;
        this->strand = strand;
        this->start = start.toUInt();
        this->end = end.toUInt();
    }
};

//! this function maps  a genomic region to genomic features
void getRegMapping(region &reg, database &anno)
{
    QVector<databaseItem*> mapping;
    // get the full recursive mapping but using priorites
    if (reg.strand == ".") {
        mapping = anno.bestRmapRange(reg.chrom, reg.start, reg.end);
    } else {
        mapping = anno.bestRmapRangeStrand(reg.chrom, reg.start, reg.end, reg.strand);
    }

    // now fill in the mappings and the best mapping variables of reg
    if (mapping.isEmpty()) {
        reg.bestmappingfeatures.clear();
    }
    else {
        foreach (databaseItem* element, mapping) {
            reg.bestmappingfeatures << element;
        }
    }
}

//! this function loads the regions and maps them to the annotation
bool loadRegions(QList<region> &regions, QString &regFile, database &anno)
{
    // variables
    bool rval = false;
    region reg;
    int regCounter = 0;


    // read the file and do all the things
    QFile file(regFile);
    if (!file.open(QFile::ReadOnly | QFile::Text)) {
        std::cerr << "Error: Cannot read file " << qPrintable(regFile)
                  << ": " << qPrintable(file.errorString())
                  << std::endl;
        rval = false;
    }
    else {
        QTextStream in(&file);
        in.setCodec("UTF-8");

        QString curline;
        QString message;
        while(!in.atEnd()) {
            curline = in.readLine().trimmed();
            QStringList fields = curline.split('\t');
            // set the original region
            reg.setInit(fields[0], fields[1], fields[2], fields[3]);
            if (fields.length() > 4) {
                for (int i = 4; i < fields.size(); ++i) {
                    reg.otherData << fields.at(i);
                }
            }
            // do the mapping
            getRegMapping(reg, anno);
            // append the nucleotide
            regions.append(reg);
            // just some time info
            ++regCounter;
            if ( (regCounter % 5000000) == 0 ) {
                message = QObject::tr("Processed %1 regions that were not skipped").arg(QString::number(regCounter));
                anno.print_time(message);
            }
        }
        message = QObject::tr("Processed %1 regions that were not skipped").arg(QString::number(regCounter));
        anno.print_time(message);

        file.close();
        if (file.error() != QFile::NoError) {
            std::cerr << "Error: Cannot read file " << qPrintable(regFile)
                      << ": " << qPrintable(file.errorString())
                      << std::endl;
            rval = false;
        }
        else { rval = true; }

    }

    return(rval);
}

//! this function writes the nucleotides to a file
bool writeRegions(QList<region> &regions, QString &fileName)
{
    bool rval = false;
    QFile file(fileName);
    if (!file.open(QFile::WriteOnly | QFile::Text)) {
        std::cerr << "Error: Cannot read file " << qPrintable(fileName)
                  << ": " << qPrintable(file.errorString())
                  << std::endl;
        rval = false;
    }
    else {
        QTextStream out(&file);
        out.setCodec("UTF-8");
        int i;

        // File writing
        foreach (region reg, regions) {
            out << reg.chrom << '\t' <<
                   reg.strand << '\t' <<
                   reg.start << '\t' <<
                   reg.end << '\t';
            if (reg.bestmappingfeatures.isEmpty()) {
                out << "none,intergenic";
            }
            else {
                out << reg.bestmappingfeatures.at(0)->topParent()->data(0).toString() << ',' << reg.bestmappingfeatures.at(0)->data(5).toString(); //the locusname and the feature
                for (i = 1; i < reg.bestmappingfeatures.count(); ++i) {
                    if (reg.bestmappingfeatures.at(i)->topParent()->data(0).toString() != reg.bestmappingfeatures.at(i-1)->topParent()->data(0).toString()) { // write only one variant
                        out << '|' << reg.bestmappingfeatures.at(i)->topParent()->data(0).toString() << ',' << reg.bestmappingfeatures.at(i)->data(5).toString(); //the locusname and the feature
                    }
                }
            }
            if (!reg.otherData.isEmpty()) {
                foreach (QString entry, reg.otherData) {
                    out << '\t' << entry;
                }
            }
            out << '\n';
        }
        out.flush();

        file.close();
        if (file.error() != QFile::NoError) {
            std::cerr << "Error: Cannot write file " << qPrintable(fileName)
                      << ": " << qPrintable(file.errorString())
                      << std::endl;
            rval = false;
        }
        else { rval = true; }
    }
    return(rval);
}



//! calculate some basic stats
void calculate_stats(QList<region> &regions)
{
    QMap<QString, float> stats; // the pair contains: context and a bestmappingfeature
    QString curkey;
    float weight;
    foreach (region reg, regions) {
        if (reg.bestmappingfeatures.isEmpty()) {
            curkey = "intergenic";
            if (!stats.contains(curkey)) { stats.insert(curkey, 0); }
            ++stats[curkey];
        }
        else {
            weight = 1/static_cast<float>(reg.bestmappingfeatures.count());
            foreach (databaseItem* element, reg.bestmappingfeatures) {
                curkey = element->data(5).toString();
                if (!stats.contains(curkey)) { stats.insert(curkey, 0); }
                stats[curkey] += weight;
            }
        }
    }
    std::cout << "feature" << '\t' << "counts" << std::endl;
    QList<QString> allKeys = stats.keys();
    foreach (curkey, allKeys) {
        std::cout << curkey.toStdString() << '\t' << stats[curkey] << std::endl;
    }
}


//! the function that is called if something with the arguments is wrong
void usage() {
    std::cerr << std::endl <<
                 "Usage is ./programname [-R regions] [-A annotation] [-O results]" << std::endl << std::endl << std::endl <<
                 "Notes" << std::endl << std::endl << '\t' <<
                 "the annotation file can be obtained via p502 format wizard (xml)." << std::endl << std::endl << '\t' <<
                 "the regions file contains a tab-separated line with the fields:" << std::endl << '\t' << '\t' <<
                 "chromosome" << std::endl << '\t' << '\t' <<
                 "strand" << std::endl << '\t' << '\t' <<
                 "start" << std::endl << '\t' << '\t' <<
                 "end" << std::endl << std::endl << '\t' <<
                 "IMPORTANT: remove any kind of header before" << std::endl << std::endl << '\t' <<
                 "IMPORTANT: coordinates need to be zero-based" << std::endl << std::endl << '\t' <<
                 "IMPORTANT: strand can be '.' which means that the mapping will be strand-UNspecific" << std::endl << std::endl << std::flush;
    exit(8);
}

int c;
extern char *optarg;
extern int optind, optopt, opterr;

int main(int argc, char *argv[])
{
    //QCoreApplication a(argc, argv);

    // check arguments
    if (argc < 4) { usage(); }

    // variables from command line
    QString regionfile = "";
    QString annofile = "";
    QString resultfile = "";
    while ((c = getopt(argc, argv, "R:A:O:")) != -1) {
        switch(c) {
        case 'R':
            regionfile = optarg;
            break;
        case 'A':
            annofile = optarg;
            break;
        case 'O':
            resultfile = optarg;
            break;
        case ':':
            std::cerr << "some stuff not specified" << std::endl << std::flush;
            break;
        case '?':
            std::cerr << "unknown argument" << std::endl << std::flush;
            usage();
        }
    }

    // create the database
    QVector<QVariant> headers;
    headers << "Sname" << "Schrom" << "Sstrand" << "Ustart" << "Uend" << "Sfeature" << "SassembledFeature" << "Upriority";
    database anno(headers);
    anno.print_time("START");
    if ( anno.readData(annofile) ) { anno.print_time("annotation loaded"); }
    else {
        std::cerr << "ERROR: could not initialize database" << std::endl << std::flush;
        exit(8);
    }

    // load the regions
    QList<region> regions;
    regions.reserve(50000000);
    if ( loadRegions(regions, regionfile, anno) ) { anno.print_time("mapping obtained"); }
    else { anno.print_time("ERROR: could not load the nucleotides"); exit(8); }
    if ( writeRegions(regions, resultfile) ) { anno.print_time("nucleotides written"); }
    else { anno.print_time("ERROR: could not write the nucleotides"); exit(8);}

    // calculate some basic stats
    calculate_stats(regions);

    // give the end
    anno.print_time("END");

    //return a.exec();
    return(0);
}





