//
// Created by Jelmer Mulder on 11/05/15.
//

#ifndef LANA_OPTIONS_H
#define LANA_OPTIONS_H

// TODO: Move this to Lana.h? Couldn't get this to work...
/// Options type
struct Options
{
    int _minCliqueSize;
    int _nMaxCliques;
    int _nMaxBlueEdges;
    bool _printProductVector;
    bool _removeAutomorphisms;
    std::string* _hist_file_name;


    Options()
            :_minCliqueSize(0)
            , _nMaxCliques(0)
            , _nMaxBlueEdges(0)
            , _printProductVector(false)
            , _removeAutomorphisms(true)
            , _hist_file_name()
    {
    }
};

#endif //LANA_OPTIONS_H
