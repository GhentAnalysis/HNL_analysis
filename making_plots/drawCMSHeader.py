
import ROOT
from ROOT import TPad, TLatex

def drawCMSHeader( pad, lumi_text, extra_text, flav_text ):

    pad.cd()

    top_margin = pad.GetTopMargin();
    header_offset = top_margin * 0.2

    left_margin = pad.GetLeftMargin()*1.27;

    header = TLatex( left_margin, 0.8, 'CMS' )
    header.SetNDC()
    header.SetTextAngle( 0 )
    header.SetTextColor( ROOT.kBlack )
    header.SetTextFont( 61 )
    header.SetTextAlign( 11 )
    cmsLabelSize = top_margin * 1.3
    header.SetTextSize( cmsLabelSize )
    cms_x_position = header.GetXsize()
    header.DrawLatex( left_margin, 0.8, 'CMS' )

    
    extra_text_size = cmsLabelSize * 0.76
    header.SetTextFont( 52 )
    header.SetTextSize( extra_text_size )
    header.DrawLatex( left_margin + 1.2 * cms_x_position, 1 - top_margin + header_offset, extra_text )
    

    lumi_text_size = top_margin * 0.7
    header.SetTextFont(42);
    header.SetTextAlign(31);
    header.SetTextSize( lumi_text_size )
    right_margin = pad.GetRightMargin();
    header.DrawLatex( 1. - right_margin, 1. - top_margin + header_offset, lumi_text )

    flav_text_size = top_margin * 0.75
    header.SetTextFont(42);
    header.SetTextAlign(31);
    header.SetTextSize( flav_text_size )
    right_margin = pad.GetRightMargin();
    header.DrawLatex( 0.7 - right_margin, 1. - top_margin + header_offset, flav_text )
